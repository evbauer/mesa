"""Bounded native HELM benchmark against frozen FreeEOS case neighborhoods."""

import argparse
import hashlib
import json
import os
import resource
import socket
import subprocess
from datetime import datetime
from pathlib import Path

import numpy as np

from native import NativeHelm, REAL, decimal_cv, response_jump, thermal_jacobian


CODE = Path(__file__).resolve().parent
WORKSPACE = CODE.parents[4]
MESA = WORKSPACE / 'mesa'
REPORT = WORKSPACE / 'freeeos_biquintic_progress'
MANIFEST = REPORT / 'experiments/038_solver_manifest/queries_v2.json'
ARCHIVE = REPORT / 'experiments/014_active_domain/active_faces.npz'
MASSFILE = MESA / 'eos/eosFreeEOS_builder/mass_frac_solar.txt'
EPSILONS = (1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def distribution(values):
    values = np.abs(np.asarray(values, dtype=float))
    return dict(zip(('count', 'min', 'median', 'p90', 'max'),
                    (len(values), float(np.min(values)), float(np.median(values)),
                     float(np.quantile(values, .9)), float(np.max(values)))))


def environment(output):
    result = dict(os.environ)
    sdk = WORKSPACE / 'mesasdk-x86_64-linux-26.6.1'
    result.update(MESA_DIR=str(MESA), MESASDK_ROOT=str(sdk),
                  MESA_CACHES_DIR=str(output / 'runtime_cache'),
                  MESA_TEMP_CACHES_DIR=str(output / 'runtime_temp_cache'),
                  LD_LIBRARY_PATH=f'{sdk}/lib:{sdk}/lib64:' + result.get('LD_LIBRARY_PATH', ''),
                  PYTHONDONTWRITEBYTECODE='1')
    result.pop('MESA_TEMP_CACHES_DISABLE', None)
    for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
        result[name] = '2'
    private_cache = output / 'runtime_cache/eosDT_cache'
    private_cache.mkdir(parents=True, exist_ok=True)
    (output / 'runtime_temp_cache').mkdir(exist_ok=True)
    installed = MESA / 'data/eosDT_data/cache/helm_table.bin'
    link = private_cache / installed.name
    if not link.exists():
        link.symlink_to(installed)
    if not link.is_symlink() or link.resolve() != installed.resolve():
        raise ValueError('private HELM cache must reuse the installed cache')
    return result


def run_kernel(output, queries, label, env):
    query_path = output / f'{label}_queries.txt'
    query_path.write_text(''.join(f'{float(density):.17e} {float(temperature):.17e}\n'
                                 for density, temperature in queries))
    command = [str(output / 'build/helm_probe'), str(MESA), str(MASSFILE), str(output)]
    with query_path.open() as source, (output / f'{label}.tsv').open('w') as target, (
            output / f'{label}.stderr').open('w') as error:
        subprocess.run(command, stdin=source, stdout=target, stderr=error,
                       cwd=output, env=env, check=True, timeout=180)
    if len(queries) == 0:
        return None
    result = np.loadtxt(output / f'{label}.tsv', skiprows=1, ndmin=2)
    if len(result) != len(queries) or not np.all(np.isfinite(result)) or np.any(result[:, 2] != 0):
        raise ValueError('invalid native HELM probe output')
    return result


def select_seeds(manifest):
    seeds = [dict(case, group='frozen_FreeEOS_24') for case in manifest['cases']]
    with np.load(ARCHIVE) as archive:
        eligible = {}
        for direction in ('u', 'tau'):
            probe = archive[f'{direction}_probe']
            for index in np.flatnonzero((probe[:, 2] == 0) & (probe[:, 3] >= 1 - 1e-12)):
                density, temperature = map(float, probe[index, :2])
                eligible.setdefault((temperature, density), (direction, int(index)))
        ranked = sorted(eligible)
    for temperature in (8.3, 8.8, 9.3):
        for density in (-2.0, 2.0, 6.0, 10.0):
            seeds.append({'case_id': f'HELM_context_{len(seeds)-23:02d}', 'group': 'HELM_context_12',
                          'logRho': density, 'logT': temperature,
                          'selection': 'predeclared Cartesian context lattice; current default-MESA routing measured separately'})
    return seeds, len(ranked)


def select_faces(table, seeds, ye):
    faces = []
    for seed in seeds:
        density = REAL(float(10**seed['logRho']))
        temperature = REAL(float(10**seed['logT']))
        din = REAL(ye) * density
        cell = table.locate(din, temperature)
        for direction, axis_number, axis in (('d', 0, table.density), ('T', 1, table.temperature)):
            for boundary in (cell[axis_number], cell[axis_number] + 1):
                if boundary <= 0 or boundary >= len(axis)-1:
                    raise ValueError('selected face reaches HELM boundary')
                point = [din, temperature]
                point[axis_number] = axis[boundary]
                left, right = list(cell), list(cell)
                left[axis_number], right[axis_number] = boundary-1, boundary
                tangential_number = 1 - axis_number
                tangential_axis = (table.density, table.temperature)[tangential_number]
                nearest = int(np.argmin(abs(np.log(point[tangential_number]/tangential_axis))))
                adjusted = abs(np.log(point[tangential_number]/tangential_axis[nearest])) < 1e-12
                if adjusted:
                    if nearest >= len(tangential_axis)-1:
                        raise ValueError('tangential knot reaches HELM boundary')
                    point[tangential_number] = (tangential_axis[nearest] + tangential_axis[nearest+1])/2
                    left[tangential_number] = nearest
                    right[tangential_number] = nearest
                faces.append({'face_id': f"{seed['case_id']}_{direction}_{boundary}",
                              'seed_id': seed['case_id'], 'group': seed['group'], 'direction': direction,
                              'din': float(point[0]), 'rho': float(point[0]/REAL(ye)), 'T': float(point[1]),
                              'logRho': float(np.log10(point[0]/REAL(ye))),
                              'logT': float(np.log10(point[1])), 'left_cell': left, 'right_cell': right,
                              'tangential_knot_moved_to_upper_linear_midpoint': bool(adjusted),
                              'displacement_dex': [float(np.log10(point[0]/REAL(ye))) - seed['logRho'],
                                                   float(np.log10(point[1])) - seed['logT']],
                              'linear_normal_widths': [float(axis[boundary]-axis[boundary-1]),
                                                       float(axis[boundary+1]-axis[boundary])]})
    return faces


def verify_source_corners(table, faces, output):
    indices = set()
    density_count = len(table.density)
    for face in faces:
        for density_index, temperature_index in (face['left_cell'], face['right_cell']):
            for density_offset in (0, 1):
                for temperature_offset in (0, 1):
                    indices.add((temperature_index+temperature_offset)*density_count + density_index+density_offset)
    largest = max(indices)
    mismatches = []
    checked = 0
    with (MESA / 'data/eosDT_data/helm_table.dat').open() as stream:
        for index, line in enumerate(stream):
            if index > largest:
                break
            if index not in indices:
                continue
            temperature_index, density_index = divmod(index, density_count)
            values = np.array([field[density_index, temperature_index] for field in table.fields])
            source = np.fromstring(line.replace('D', 'E'), sep=' ')
            if len(source) != 9 or not np.array_equal(source, values):
                mismatches.append(index)
            checked += 1
    result = {'unique_corners': checked, 'fields_per_corner': 9,
              'cache_matches_printed_source_bitwise': not mismatches,
              'mismatched_flat_indices': mismatches,
              'scope': 'all nine potential slots at every corner used in this benchmark'}
    dump(output / 'source_corner_check.json', result)
    if checked != len(indices) or mismatches:
        raise ValueError('queried source/cache potential corners disagree')
    return result


def measure_faces(table, faces, ye, output, env):
    queries, query_map = [], []
    for face_index, face in enumerate(faces):
        queries.append((face['rho'], face['T']))
        query_map.append((face_index, 'center', -1, 0))
        for epsilon in EPSILONS:
            for side, sign in (('left', -1), ('right', 1)):
                point = [REAL(face['din']), REAL(face['T'])]
                axis_number = 0 if face['direction'] == 'd' else 1
                width = face['linear_normal_widths'][0 if side == 'left' else 1]
                point[axis_number] += REAL(sign * epsilon * width)
                queries.append((float(point[0]/REAL(ye)), float(point[1])))
                query_map.append((face_index, side, EPSILONS.index(epsilon), epsilon))
    kernel = run_kernel(output, queries, 'native', env)
    checks = []
    for face in faces:
        face['near_checks'] = []
        point = [REAL(face['din']), REAL(face['T'])]
        for side in ('left', 'right'):
            reference = decimal_cv(table, *point, ye, face[side+'_cell'])
            face['electron_'+side] = [float(value) for value in reference]
            face['electron_'+side+'_decimal'] = [str(value) for value in reference]
        face['float64_endpoint_error_over_Cv'] = max(float(np.max(np.abs(
            table.cv(*point, ye, face[side+'_cell'], dtype=np.float64) - table.cv(*point, ye, face[side+'_cell']))))
            for side in ('left', 'right')) / max(1, abs(face['electron_left'][0]))
        face['electron_jump'] = response_jump(face['electron_left'], face['electron_right'])
    for row, (face_index, side, epsilon_index, epsilon) in zip(kernel, query_map):
        face = faces[face_index]
        if side == 'center':
            additions = row[6:9] + row[9:12] + row[12:15]
            face['center_components'] = {name: row[start:start+3].tolist()
                                         for name, start in (('electron_kernel', 3), ('ion', 6), ('radiation', 9),
                                                             ('coulomb_by_subtraction', 12), ('total_kernel', 15))}
            face['Gamma'] = float(row[18])
            for endpoint in ('left', 'right'):
                face['total_'+endpoint] = (np.asarray(face['electron_'+endpoint]) + additions).tolist()
            face['total_jump'] = response_jump(face['total_left'], face['total_right'])
            face['thermal_jacobian'] = [thermal_jacobian(face['total_left'], face['total_right'], increment)
                                        for increment in (0, 1e-4, 1e-3, 1e-2, -1e-4, -1e-3, -1e-2)]
            continue
        actual_din = REAL(ye) * REAL(row[0])
        actual_temperature = REAL(row[1])
        independent = table.cv(actual_din, actual_temperature, ye, face[side+'_cell'])
        exact_electron = np.asarray(face['electron_'+side], dtype=REAL)
        exact_total = np.asarray(face['total_'+side], dtype=REAL)
        response_scale = max(1, abs(exact_total[0]), abs(exact_total[1]), abs(exact_total[2]))
        check = {'side': side, 'epsilon_fraction_of_linear_cell': epsilon,
                 'kernel_vs_independent_over_scale': float(np.max(np.abs(row[3:6] - independent))/response_scale),
                 'independent_electron_limit_error_over_scale': float(np.max(np.abs(independent-exact_electron))/response_scale),
                 'kernel_total_limit_error_over_scale': float(np.max(np.abs(row[15:18]-exact_total))/response_scale),
                 'Gamma': float(row[18])}
        face['near_checks'].append(check)
        checks.append(check)
    return checks


def route_states(faces, seeds, output, env):
    route_dir = output / 'routing'
    route_dir.mkdir(exist_ok=True)
    states = seeds + faces
    query_path = route_dir / 'queries.txt'
    query_path.write_text(''.join(f"{state['logRho']:.17e} {state['logT']:.17e}\n" for state in states))
    command = ['python3', str(CODE.parent/'mesa_probe/run_probe.py'), '--output', str(route_dir)]
    with query_path.open() as source, (route_dir/'probe.tsv').open('w') as target, (
            route_dir/'probe.stderr').open('w') as error:
        subprocess.run(command, stdin=source, stdout=target, stderr=error, env=env,
                       cwd=route_dir, check=True, timeout=180)
    rows = np.loadtxt(route_dir/'probe.tsv', skiprows=1, ndmin=2)
    if len(rows) != len(states) or np.any(rows[:, 2] != 0) or not np.isfinite(rows).all():
        raise ValueError('default-MESA route query failed')
    for state, row in zip(states, rows):
        state['default_MESA'] = {'frac_HELM': float(row[3]), 'frac_FreeEOS': float(row[5]),
                                 'all_fractions': dict(zip(('HELM', 'OPAL_SCVH', 'FreeEOS', 'PC', 'Skye', 'CMS', 'ideal', 'WDM_AV'), row[3:11].tolist())),
                                 'Cv': float(row[11]), 'Cv_x': float(row[12]), 'Cv_t': float(row[13])}
    return len(rows)


def matched_seed_states(seeds, output, env):
    rows = run_kernel(output, [(10**seed['logRho'], 10**seed['logT']) for seed in seeds], 'seeds', env)
    for seed, row in zip(seeds, rows):
        seed['HELM_at_exact_seed'] = dict(zip(('Cv', 'Cv_x', 'Cv_t'), row[15:18].tolist()))
        if seed['group'] == 'frozen_FreeEOS_24':
            left = [seed['split_left'][key] for key in ('Cv', 'Cv_x', 'Cv_t')]
            right = [seed['split_right'][key] for key in ('Cv', 'Cv_x', 'Cv_t')]
            seed['FreeEOS_jump'] = response_jump(left, right)
            seed['FreeEOS_thermal_jacobian'] = [thermal_jacobian(left, right, increment)
                                               for increment in (0, 1e-4, 1e-3, 1e-2, -1e-4, -1e-3, -1e-2)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=REPORT/'experiments/048_helm_control')
    args = parser.parse_args()
    if socket.gethostname().split('.')[0] != 'rzwhippet14' or os.environ.get('SLURM_JOB_ID') != '5528075':
        raise RuntimeError('TASK048 placement/allocation guard failed')
    resource.setrlimit(resource.RLIMIT_AS, (8*1024**3, 8*1024**3))
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    initial_paths = {}
    for line in (output/'initial_sha256.txt').read_text().splitlines():
        digest, relative = line.split(maxsplit=1)
        initial_paths[relative] = digest
    print('Checking pinned dependencies', flush=True)
    for relative, digest in initial_paths.items():
        if sha256(WORKSPACE/relative) != digest:
            raise RuntimeError('initial dependency changed: ' + relative)
    if sha256(MANIFEST) != 'afc6fe676ba13bf858fc0a611268402b3663f4ad1e0ff991049713fd78788f19':
        raise ValueError('frozen manifest changed')
    if sha256(MASSFILE) != 'bd62922804e073a2a410480994b8a28e85b0bc53e469729d6ed5136ea9da10a6':
        raise ValueError('solar composition pattern changed')
    env = environment(output)
    run_kernel(output, [], 'initialize', env)
    composition = dict(zip(('X', 'Y', 'Z', 'abar', 'zbar', 'Ye', 'crad', 'kerg', 'avo'),
                           np.loadtxt(output/'composition.tsv', skiprows=1).tolist()))
    table = NativeHelm(MESA/'data/eosDT_data/cache/helm_table.bin',
                       output/'density_axis.txt', output/'temperature_axis.txt')
    manifest = json.loads(MANIFEST.read_text())
    if sha256(ARCHIVE) != manifest['provenance']['archive_sha256']:
        raise ValueError('frozen 014 source archive changed')
    seeds, eligible_count = select_seeds(manifest)
    faces = select_faces(table, seeds, composition['Ye'])
    dump(output/'sampling.json', {'seeds': seeds, 'faces': faces,
                                  'archived_pure_HELM_unique_eligible': eligible_count})
    print(f'Frozen sampling: {len(seeds)} seeds, {len(faces)} nearby native faces', flush=True)
    source_check = verify_source_corners(table, faces, output)
    checks = measure_faces(table, faces, composition['Ye'], output, env)
    matched_seed_states(seeds, output, env)
    routing_count = route_states(faces, seeds, output, env)
    print('Native traces, convergent offsets, and unchanged default-MESA routing complete', flush=True)
    summaries = {}
    for group in ('frozen_FreeEOS_24', 'HELM_context_12'):
        for direction in ('d', 'T'):
            selected = [face for face in faces if face['group'] == group and face['direction'] == direction]
            summaries[group+'_'+direction] = {
                component+'_'+field: distribution([face[component+'_jump'][field]['signed_scaled'] for face in selected])
                for component in ('electron', 'total') for field in ('Cv_x', 'Cv_t')}
            summaries[group+'_'+direction]['routing'] = {
                'pure_HELM': sum(face['default_MESA']['frac_HELM'] >= 1-1e-12 for face in selected),
                'pure_FreeEOS': sum(face['default_MESA']['frac_FreeEOS'] >= 1-1e-12 for face in selected)}
            summaries[group+'_'+direction]['thermal_delta_1e-3'] = distribution([
                face['thermal_jacobian'][2]['signed_scaled'] for face in selected])
    freeeos = [seed for seed in seeds if seed['group'] == 'frozen_FreeEOS_24']
    summaries['FreeEOS_frozen_24'] = {field: distribution([
        seed['FreeEOS_jump'][field]['signed_scaled'] for seed in freeeos]) for field in ('Cv_x', 'Cv_t')}
    summaries['FreeEOS_frozen_24']['thermal_delta_1e-3'] = distribution([
        seed['FreeEOS_thermal_jacobian'][2]['signed_scaled'] for seed in freeeos])
    validation = {'source_corners': source_check, 'routing_queries': routing_count,
                  'kernel_queries': len(faces)*(1 + 2*len(EPSILONS)) + len(seeds),
                  'kernel_vs_extended_precision': distribution([check['kernel_vs_independent_over_scale'] for check in checks]),
                  'near_limit_by_epsilon': {str(epsilon): distribution([
                      check['kernel_total_limit_error_over_scale'] for check in checks
                      if check['epsilon_fraction_of_linear_cell'] == epsilon]) for epsilon in EPSILONS},
                  'longdouble_mantissa_bits': np.finfo(REAL).nmant,
                  'min_distance_Gamma_to_one': min(abs(face['Gamma'] - 1) for face in faces),
                  'Coulomb_branch_unchanged_at_all_offsets': all((check['Gamma'] >= 1) == (face['Gamma'] >= 1)
                      for face in faces for check in face['near_checks']),
                  'all_initial_dependencies_unchanged': all(sha256(WORKSPACE/path) == digest
                                                           for path, digest in initial_paths.items())}
    result = {'schema': 'NATIVE_HELM_CONTROL_V1', 'created': datetime.now().isoformat(),
              'host': socket.gethostname(), 'job': os.environ['SLURM_JOB_ID'],
              'provenance': {'initial_dependencies_sha256': initial_paths,
                  'archive_sha256': sha256(ARCHIVE),
                  'source_code_sha256': {name: sha256(CODE/name) for name in ('native.py', 'run_control.py')},
                  'private_driver_sha256': sha256(CODE/'helm_probe.f90')},
              'composition': composition, 'sampling': {'frozen_FreeEOS_seeds': 24, 'HELM_context_seeds': 12,
                  'native_faces_per_seed': 4, 'method': 'two bracketing native T faces and two bracketing native d faces; seed tangent retained except near-native knots (<1e-12 ln) move to the upper linear midpoint',
                  'HELM_context_selection': 'Cartesian logT=(8.3,8.8,9.3), logRho=(-2,2,6,10), fixed before jump/routing results',
                  'eligible_HELM_states': eligible_count},
              'coordinates': {'axes': '100 intervals/dex; polynomial coordinates linear T and d=Ye*rho',
                  'Cv': '-Ye*T*f_TT', 'Cv_x': '-Ye*T*d*f_dTT', 'Cv_t': '-Ye*T*(f_TT+T*f_TTT)'},
              'conventions': {'radiation': 'MESA crad, fully restored Planck radiation',
                  'ions': 'fully ionized MESA abar/zbar; native ideal ions',
                  'Coulomb': 'native corrections enabled; both cutoffs -1e99 as installed defaults; derivatives from total minus ep/ion/rad',
                  'jump_sign': 'left-minus-right; summaries give absolute values',
                  'response_scale': 'max(1,abs(Cv_left),abs(Cv_right),abs(partial_left),abs(partial_right))',
                  'consumer': 'constant-density eps=-T*Cv*delta_lnT/dt; dt=1; d(eps)/dlnT=-T*((1+delta_lnT)*Cv+delta_lnT*Cv_t)',
                  'consumer_scale': 'max(T*abs(Cv_left/right),abs(J_left/right)); values stored divided by common T',
                  'total_endpoints': '80-digit electron endpoint plus common exact-face native ion/rad/Coulomb additions; checked against kernel limits'},
              'validation': validation, 'summaries': summaries, 'seeds': seeds, 'faces': faces}
    dump(output/'results.json', result)
    print(json.dumps({'summaries': summaries, 'validation': validation}, indent=2), flush=True)
    if not validation['all_initial_dependencies_unchanged']:
        raise RuntimeError('protected dependency changed during benchmark')


if __name__ == '__main__':
    main()
