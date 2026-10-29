"""Table-implied pressure-strip audit, not a verified FreeEOS source constant."""

import argparse
import csv
import hashlib
import json
import platform
import sys
from pathlib import Path

import numpy as np


WORKSPACE = Path(__file__).resolve().parents[5]
TABLE = WORKSPACE / 'mesa/data/eosFreeEOS_data/mesa-FreeEOS_02z70x.data'
OUTPUT = WORKSPACE / 'freeeos_biquintic_progress/experiments/005_radiation_audit'
EXPECTED_SHA = 'e6ca7ce27bee8197bb257a8c4886a06366e6951e58958fbed3735c877ee68502'
CRAD = 8.0 * np.pi**5 * (1.380649e-16)**4 / (
    15.0 * (6.62607015e-27)**3 * (2.99792458e10)**3)
LN10 = np.log(10.0)


def print_half_unit(token):
    mantissa, exponent = token.lower().split('e')
    decimal_places = len(mantissa.split('.')[1])
    return 0.5 * 10.0**(int(exponent) - decimal_places)


def read_table(path):
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != EXPECTED_SHA:
        raise ValueError(f'Unexpected input SHA-256: {digest}')
    records, bounds, block_labels = [], [], []
    pending_block = False
    with path.open() as stream:
        next(stream)
        header = np.fromstring(next(stream), sep=' ')
        for record in stream:
            if 'logQ = logRho' in record:
                pending_block = True
                continue
            if pending_block and record.strip():
                block_labels.append(float(record))
                pending_block = False
                continue
            tokens = record.split()
            if len(tokens) == 22 and tokens[0][0] in '+-0123456789':
                records.append([float(token) for token in tokens])
                bounds.append([print_half_unit(token) for token in tokens])
    if header.size != 11 or header[0] != 51:
        raise ValueError('Expected format-51 header')
    shape = (int(header[7]), int(header[3]), 22)
    if len(records) != shape[0] * shape[1] or len(block_labels) != shape[0]:
        raise ValueError('Incomplete table')
    return header, np.asarray(records).reshape(shape), np.asarray(bounds).reshape(shape), digest


def infer_fields(rows, rounding):
    temperature = 10.0**rows[..., 0]
    rho = 10.0**rows[..., 18]
    pgas = 10.0**rows[..., 1]
    energy = 10.0**rows[..., 2]
    entropy = 10.0**rows[..., 3]
    chirho, chit = rows[..., 4], rows[..., 5]
    cp, cv = rows[..., 6], rows[..., 7]
    erho, st, srho = rows[..., 8], rows[..., 9], rows[..., 10]
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        prad = CRAD * temperature**4 / 3.0
        current = pgas + prad
        maxwell = -rho**2 * temperature * srho / chit
        energy_balance = rho**2 * erho / (1.0 - chit)
        strip_term = 3.0 * pgas / temperature**4
        maxwell_term = 3.0 * maxwell / temperature**4
        energy_term = 3.0 * energy_balance / temperature**4
        a_maxwell = maxwell_term - strip_term
        a_energy = energy_term - strip_term
        rho_error = LN10 * rounding[..., 18]
        temperature_error = LN10 * rounding[..., 0]
        pgas_error = LN10 * rounding[..., 1]
        maxwell_round = abs(maxwell_term) * (
            2 * rho_error + 3 * temperature_error
            + rounding[..., 10] / abs(srho) + rounding[..., 5] / abs(chit))
        energy_round = abs(energy_term) * (
            2 * rho_error + 4 * temperature_error
            + rounding[..., 8] / abs(erho) + rounding[..., 5] / abs(1 - chit))
        strip_round = abs(strip_term) * (pgas_error + 4 * temperature_error)
        arithmetic_maxwell = 20 * np.finfo(float).eps * (abs(maxwell_term) + abs(strip_term))
        arithmetic_energy = 20 * np.finfo(float).eps * (abs(energy_term) + abs(strip_term))
        maxwell_round += strip_round + arithmetic_maxwell
        energy_round += strip_round + arithmetic_energy
        erho_proxy = -rho**2 * erho / temperature**4
        erho_proxy_round = abs(erho_proxy) * (
            2 * rho_error + 4 * temperature_error + rounding[..., 8] / abs(erho)
            + 20 * np.finfo(float).eps)
        linear_density_pressure_proxy = 3 * current * (1 - chirho) / temperature**4
        results = {
            'rho': rho, 'T': temperature, 'Pgas': pgas, 'Prad_over_Pgas_current': prad / pgas,
            'P_current': current, 'P_maxwell': maxwell, 'P_energy': energy_balance,
            'a_maxwell': a_maxwell, 'a_energy': a_energy,
            'a_maxwell_over_crad_minus1': a_maxwell / CRAD - 1,
            'a_energy_over_crad_minus1': a_energy / CRAD - 1,
            'a_route_difference_over_crad': (a_maxwell - a_energy) / CRAD,
            'a_maxwell_rounding_bound_over_crad': maxwell_round / CRAD,
            'a_energy_rounding_bound_over_crad': energy_round / CRAD,
            'a_route_difference_in_rounding_bounds': abs(a_maxwell - a_energy) / (maxwell_round + energy_round),
            'a_maxwell_current_difference_in_rounding_bounds': abs(a_maxwell - CRAD) / maxwell_round,
            'a_energy_current_difference_in_rounding_bounds': abs(a_energy - CRAD) / energy_round,
            'a_erho_proxy_over_crad_minus1': erho_proxy / CRAD - 1,
            'a_erho_proxy_rounding_bound_over_crad': erho_proxy_round / CRAD,
            'a_linear_density_pressure_proxy_over_crad_minus1': linear_density_pressure_proxy / CRAD - 1,
            'a_erho_minus_linear_density_pressure_proxy_over_crad': (erho_proxy - linear_density_pressure_proxy) / CRAD,
            'erho_proxy_excess_pressure_over_printed_Pgas': (erho_proxy / CRAD - 1) * prad / pgas,
            'a_srho_proxy_over_crad_minus1': -0.75 * rho**2 * srho / (temperature**3 * CRAD) - 1,
            'a_cv_proxy_over_crad_minus1': rho * cv / (4 * temperature**3 * CRAD) - 1,
            'a_e_proxy_over_crad_minus1': rho * energy / (temperature**4 * CRAD) - 1,
            'a_s_proxy_over_crad_minus1': 0.75 * rho * entropy / (temperature**3 * CRAD) - 1,
            'maxwell_identity_current': maxwell / current - 1,
            'energy_identity_current': rho**2 * erho / current + chit - 1,
            'cv_entropy_identity': temperature * st / cv - 1,
            'chiT': chit, 'chiRho': chirho, 'one_minus_chiT_abs': abs(1 - chit),
            'a_maxwell_subtraction_condition': (abs(maxwell_term) + abs(strip_term)) / abs(a_maxwell),
            'a_energy_subtraction_condition': (abs(energy_term) + abs(strip_term)) / abs(a_energy),
            'total_over_gas_pressure_condition': current / pgas,
            'inverse_chiRho_condition': 1 / chirho,
            'energy_chiT_absolute_condition': 1 / abs(1 - chit),
        }
        for name, pressure in [('current', current), ('maxwell', maxwell), ('energy', energy_balance)]:
            expected_cp = cv + pressure * chit**2 / (rho * temperature * chirho)
            results[f'cp_identity_{name}'] = expected_cp / cp - 1
            results[f'P_chiRho_over_Pgas_{name}'] = pressure * chirho / pgas
    return results


def clean_mask(rows, header, radius):
    logq = header[8] + np.arange(rows.shape[0]) * header[10]
    expected_t = header[4] + np.arange(rows.shape[1]) * header[6]
    coordinate_error = rows[..., 18] - (logq[:, None] + 2 * rows[..., 0] - 12)
    valid = (rows[..., 17] == 0) & (rows[..., 0] >= 3.2)
    valid &= np.all(np.isfinite(rows), axis=-1)
    valid &= np.all(rows[..., [4, 5, 6, 7]] > 0, axis=-1)
    valid &= abs(coordinate_error) < 1e-10
    valid &= abs(rows[..., 0] - expected_t[None, :]) < 1e-10
    valid[520, :] = False
    clean = np.zeros_like(valid)
    if radius == 0:
        return valid, valid.copy(), coordinate_error
    interior = np.ones((valid.shape[0] - 2 * radius, valid.shape[1] - 2 * radius), dtype=bool)
    for density_offset in range(2 * radius + 1):
        for temperature_offset in range(2 * radius + 1):
            interior &= valid[density_offset:density_offset + interior.shape[0],
                              temperature_offset:temperature_offset + interior.shape[1]]
    clean[radius:-radius, radius:-radius] = interior
    return valid, clean, coordinate_error


def distribution(values):
    values = np.asarray(values)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return {'count': int(values.size), 'finite': 0}
    quantiles = np.quantile(finite, [0, .01, .5, .99, 1])
    absolute = np.quantile(abs(finite), [.5, .99, 1])
    return {'count': int(values.size), 'finite': int(finite.size),
            **dict(zip(['min', 'p01', 'median', 'p99', 'max'], map(float, quantiles))),
            **dict(zip(['abs_median', 'abs_p99', 'abs_max'], map(float, absolute)))}


def summarize(fields, mask):
    return {name: distribution(values[mask]) for name, values in fields.items()}


def write_csv(path, columns):
    with path.open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(columns)
        for record in zip(*columns.values()):
            writer.writerow([format(float(value), '.17g') for value in record])


def run(args):
    header, rows, rounding, digest = read_table(args.table)
    valid, clean, coordinate_error = clean_mask(rows, header, args.stencil_radius)
    fields = infer_fields(rows, rounding)
    ratio = fields['Prad_over_Pgas_current']
    selected = clean & (ratio >= args.min_ratio)
    selected &= fields['one_minus_chiT_abs'] >= args.min_one_minus_chit
    selected &= (fields['P_maxwell'] > 0) & (fields['P_energy'] > 0)
    selected &= np.isfinite(fields['a_maxwell']) & np.isfinite(fields['a_energy'])
    if not selected.any():
        raise ValueError('No conditioned radiation-dominated clean nodes')
    strata = {}
    for threshold in [1e2, 1e4, 1e5, 2e5, 1e6]:
        strata[f'ratio_ge_{threshold:g}'] = summarize(fields, selected & (ratio >= threshold))
    tight_mask = selected & (rows[..., 0] >= 7.5) & (ratio >= 2e5)
    tight_stats = summarize(fields, tight_mask)
    trends = []
    for axis, coordinates, edges in [
            ('logT', rows[..., 0], np.arange(3.0, 8.51, .5)),
            ('logRho', rows[..., 18], np.arange(-20.0, 12.1, 2.0)),
            ('log10_ratio', np.log10(ratio), np.arange(2.0, 19.1, 1.0))]:
        for lower, upper in zip(edges[:-1], edges[1:]):
            mask = selected & (coordinates >= lower) & (coordinates < upper)
            if mask.any():
                trends.append({'axis': axis, 'lower': float(lower), 'upper': float(upper),
                               'statistics': summarize(fields, mask)})
    summary = {
        'claim': 'Table-implied pressure-strip reconstruction audit; NOT a verified source constant or source accuracy.',
        'table': str(args.table.resolve()), 'table_sha256': digest, 'current_crad': float(CRAD),
        'python': platform.python_version(), 'numpy': np.__version__, 'command': sys.argv,
        'selection': {'min_logT': 3.2, 'excluded_iq_zero_based': [520],
                      'coordinate_tolerance_dex': 1e-10, 'coordinates': 'printed logRho and logT',
                      'native_stencil_radius': args.stencil_radius, 'min_Prad_over_Pgas_current': args.min_ratio,
                      'min_abs_one_minus_chiT': args.min_one_minus_chit,
                      'total_nodes': int(valid.size), 'valid_native_nodes': int(valid.sum()),
                      'clean_stencil_nodes': int(clean.sum()), 'selected_nodes': int(selected.sum())},
        'rounding': 'Half last printed decimal unit propagated conservatively to a, plus 20*float64 epsilon per route; not solver error bars.',
        'statistics': summarize(fields, selected), 'dominance_strata': strata, 'trends': trends,
        'high_T_high_dominance_diagnostic': {
            'min_logT': 7.5, 'min_Prad_over_Pgas_current': 2e5, 'statistics': tight_stats,
            'interpretation': 'Not a fit: density-derivative energy proxy and pressure proxy under a linear-density material-pressure hypothesis; material corrections are not silently removed.'},
        'rounding_bound_exceedance_counts': {
            'current_a_maxwell': int((selected & (fields['a_maxwell_current_difference_in_rounding_bounds'] > 1)).sum()),
            'current_a_energy': int((selected & (fields['a_energy_current_difference_in_rounding_bounds'] > 1)).sum()),
            'between_routes': int((selected & (fields['a_route_difference_in_rounding_bounds'] > 1)).sum())},
        'limitations': [
            'Both identities share chiT, density and originating table; route agreement is closure, not independent accuracy.',
            'Inferred a subtracts printed Pgas from identity-implied total P; it diagnoses pressure stripping, not automatically the internal energy radiation constant.',
            'erho/srho/Cv/E/S proxies include material thermodynamics or reference gauges; high-dominance trends alone do not verify source conventions.',
            'The linear-density pressure proxy additionally assumes material pressure is locally linear in density; it is a diagnostic, not a general EOS identity.',
            'Current-CRAD ratio defines selection; no constant or other-physics fit is performed.',
            'Small chiRho can amplify material remainder and heat-capacity sensitivity even when total-P closure is good.',
        ],
    }
    OUTPUT.mkdir(parents=True, exist_ok=True)
    products = [OUTPUT / name for name in ['results.json', 'nodes.csv', 'trends.csv', 'SUMMARY.md']]
    if not args.overwrite and any(path.exists() for path in products):
        raise ValueError('Audit products already exist; use --overwrite to rerun this experiment')
    products[0].write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
    density_index, temperature_index = np.nonzero(selected)
    write_csv(products[1], {'iq': density_index, 'it': temperature_index,
                           'logQ_header': (header[8] + density_index * header[10]),
                           'logT_printed': rows[..., 0][selected], 'logRho_printed': rows[..., 18][selected],
                           'coordinate_error_dex': coordinate_error[selected],
                           **{name: values[selected] for name, values in fields.items()}})
    with products[2].open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['axis', 'lower', 'upper', 'count', 'maxwell_bias_median', 'energy_bias_median',
                         'maxwell_abs_p99', 'energy_abs_p99', 'route_abs_p99', 'erho_proxy_bias_median',
                         'srho_proxy_bias_median', 'cp_current_abs_p99'])
        for trend in trends:
            stats = trend['statistics']
            writer.writerow([trend['axis'], trend['lower'], trend['upper'], stats['a_maxwell']['count'],
                             stats['a_maxwell_over_crad_minus1']['median'], stats['a_energy_over_crad_minus1']['median'],
                             stats['a_maxwell_over_crad_minus1']['abs_p99'], stats['a_energy_over_crad_minus1']['abs_p99'],
                             stats['a_route_difference_over_crad']['abs_p99'],
                             stats['a_erho_proxy_over_crad_minus1']['median'],
                             stats['a_srho_proxy_over_crad_minus1']['median'], stats['cp_identity_current']['abs_p99']])
    lines = ['# Table-implied radiation audit', '', summary['claim'], '',
             f'Input SHA-256: `{digest}`; current CRAD: `{CRAD:.17e}`.',
             f'Selected {selected.sum()} clean native nodes, ratio >= {args.min_ratio:g}, '
             f'|1-chiT| >= {args.min_one_minus_chit:g}, stencil radius {args.stencil_radius}.', '',
             '| Minimum current Prad/Pgas | Nodes | Maxwell median a/CRAD-1 | Energy median | Route abs p99 | Cp current abs p99 |',
             '| ---: | ---: | ---: | ---: | ---: | ---: |']
    for threshold, stats in strata.items():
        count = stats['a_maxwell']['finite']
        if count:
            lines.append(f"| {threshold.removeprefix('ratio_ge_')} | {count} | "
                         f"{stats['a_maxwell_over_crad_minus1']['median']:.6e} | "
                         f"{stats['a_energy_over_crad_minus1']['median']:.6e} | "
                         f"{stats['a_route_difference_over_crad']['abs_p99']:.6e} | "
                         f"{stats['cp_identity_current']['abs_p99']:.6e} |")
    lines.extend(['', 'Rounding and rho/T trends: `results.json`, `trends.csv`; selected raw audit fields: `nodes.csv`.', '',
                  'Energy/entropy radiation proxies include material contributions; no source constant or gauge is verified.',
                  'Pressure-strip and internal-energy interpretations must be distinguished. No constant fit or other-physics removal is performed.'])
    if tight_stats['a_maxwell']['finite']:
        erho_bias = tight_stats['a_erho_proxy_over_crad_minus1']['median']
        linear_bias = tight_stats['a_linear_density_pressure_proxy_over_crad_minus1']['median']
        lines.extend(['', '## High-temperature, high-dominance diagnostic', '',
                      f"{tight_stats['a_maxwell']['finite']} nodes at logT>=7.5 and ratio>=2e5:",
                      f'- Median energy-density radiation proxy/CRAD-1: {erho_bias:.12e}.',
                      f'- Corresponding dimensional proxy: {CRAD * (1 + erho_bias):.17e}.',
                      f'- Median linear-density pressure proxy/CRAD-1: {linear_bias:.12e}.',
                      f"- Proxy difference abs p99: {tight_stats['a_erho_minus_linear_density_pressure_proxy_over_crad']['abs_p99']:.6e} CRAD.",
                      f"- Energy-proxy rounding bound p99: {tight_stats['a_erho_proxy_rounding_bound_over_crad']['p99']:.6e} CRAD.",
                      f"- Proxy excess pressure / printed Pgas median: {tight_stats['erho_proxy_excess_pressure_over_printed_Pgas']['median']:.6e}.",
                      '', 'The pressure-strip a agrees with current CRAD, but energy-density and linear-density pressure proxies',
                      'show a common positive residual far above printing estimates. This supports a distinct radiation-like',
                      'contribution in archived energy/material remainders, not a verified upstream radiation constant.',
                      'Keep current CRAD for printed-pressure reconstruction; do not replace it with the energy proxy.',
                      'Current-CRAD subtraction is therefore not proven to isolate pure material thermodynamics.'])
    products[3].write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--table', type=Path, default=TABLE)
    parser.add_argument('--min-ratio', type=float, default=100.0)
    parser.add_argument('--min-one-minus-chit', type=float, default=.25)
    parser.add_argument('--stencil-radius', type=int, default=2)
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()
    if args.min_ratio < 100 or args.min_one_minus_chit <= 0 or not 0 <= args.stencil_radius <= 10:
        parser.error('Require ratio >=100, positive denominator threshold, and stencil radius 0..10')
    run(args)


if __name__ == '__main__':
    main()
