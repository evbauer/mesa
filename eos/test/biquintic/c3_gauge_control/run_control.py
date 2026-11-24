"""Five coarse fits exposing known-gauge sensitivity and analytic-reference removal."""

import argparse
import hashlib
import json
import resource
import sys
import time
from datetime import datetime
from importlib.metadata import version

import gauge
import numpy as np
from gauge import (
    BIQUINTIC, HERE, LOW_ORDERS, analytic_gauge_roughness, recovered_states,
    remove_input_gauge, shifted_low, solve_regularized,
)
from regularized import ARCHIVED, FIELDS, LN10, TABLE, TABLE_SHA, WORKSPACE, positive_status, summary
from freeeos_potential import read_table, source_constraints

OUTPUT = WORKSPACE / 'freeeos_biquintic_progress/experiments/034_c3_gauge'
SOFT = WORKSPACE / 'freeeos_biquintic_progress/experiments/031_regularized_c3'
HARD = WORKSPACE / 'freeeos_biquintic_progress/experiments/033_first_constraint'
FROZEN_SHA = 'ba30d62f864e094bb52ad0b83f87554b8fff199c4146f1c1e1a6fc297b35edfd'
RESPONSE_PARTIALS = tuple(field + suffix for field in ['Cv','Cp','chiRho','chiT'] for suffix in ['_x','_t'])
OUTPUT_FIELDS = FIELDS + RESPONSE_PARTIALS


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def departures(state, reference):
    result = {}
    for field in OUTPUT_FIELDS:
        if field in RESPONSE_PARTIALS:
            base = field[:-2]
            scale = np.maximum(1., np.maximum(abs(reference[base]), abs(reference[field])))
        else:
            scale = abs(reference[field])
        result[field] = summary((state[field] - reference[field]) / scale)
    return result


def run(args):
    started = time.perf_counter()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    if (OUTPUT / 'results.json').exists() and not args.overwrite:
        raise ValueError('Products exist; use --overwrite')
    helper = BIQUINTIC / 'c3_spline_control/regularized.py'
    if digest(helper) != FROZEN_SHA or digest(TABLE) != TABLE_SHA:
        raise ValueError('Frozen API or pinned table mismatch')
    patch = json.loads((ARCHIVED / 'results.json').read_text())['patches']['active_chiT']
    soft_archive = np.load(SOFT / 'paired_fits_and_queries.npz')
    hard_archive = np.load(HARD / 'paired_controls.npz')
    _, logq, logt, rows = read_table(TABLE)
    low, valid, _ = source_constraints(rows, logq, logt)
    start_q, start_t, count_q, count_t = patch['selected_indices']
    grid_q, grid_t = logq[start_q:start_q+count_q:2], logt[start_t:start_t+count_t:2]
    coarse_low = low[:, :, start_q:start_q+count_q:2, start_t:start_t+count_t:2]
    if not valid[start_q:start_q+count_q:2,start_t:start_t+count_t:2].all():
        raise ValueError('Coarse source patch no longer all native valid')
    sample_q, sample_t = soft_archive['withheld_logQ'], soft_archive['withheld_logT']
    np.testing.assert_array_equal(sample_q, hard_archive['withheld_logQ'])
    np.testing.assert_array_equal(sample_t, hard_archive['withheld_logT'])
    node_q, node_t = np.meshgrid(grid_q, grid_t, indexing='ij')
    node_q, node_t = node_q.ravel(), node_t.ravel()
    coarse_rows = rows[start_q:start_q+count_q:2,start_t:start_t+count_t:2]
    a = float(np.median(10.0**coarse_rows[:, :, 2]))
    b = float(np.median(10.0**coarse_rows[:, :, 3]))
    source_files = [HERE/'gauge.py', HERE/'run_control.py', HERE/'test_gauge.py', helper,
                    BIQUINTIC/'c3_spline_control/spline_control.py', BIQUINTIC/'freeeos_potential.py']
    result = {'claim': 'Known artificial gauge sensitivity of a table-only roughness selector; no inference about original FreeEOS energy/entropy gauge, fresh-source accuracy or runtime.',
        'table_sha256': TABLE_SHA, 'frozen031_helper_sha256': FROZEN_SHA,
        'source_hashes': {str(path.relative_to(WORKSPACE)):digest(path) for path in source_files},
        'archived_hashes': {str(path.relative_to(WORKSPACE)):digest(path) for path in [ARCHIVED/'results.json',
            SOFT/'results.json', SOFT/'paired_fits_and_queries.npz', HARD/'results.json', HARD/'paired_controls.npz']},
        'versions': {'python':sys.version, 'numpy':np.__version__, 'scipy':version('scipy')},
        'command':sys.argv, 'timestamp':datetime.now().astimezone().isoformat(),
        'patch':{'label':'active_chiT','bounds_from015':patch['actual_bounds'],'coarse_shape':list(coarse_rows.shape[:2]),
                 'unknowns':1064,'native_source_nodes':len(node_q),'same_saved_odd_odd_holdouts':len(sample_q)},
        'gauges':{'a_erg_per_g':a,'b_erg_per_g_K':b,'definitions':'a=median printed retained E; b=median printed retained S; shifts a+b*T and a+100*b*T. Artificial known offsets, not original gauge estimates.',
                  'known_physical_effect':'E -> E+a, S -> S-b (or -100b); P,Cv,Cp,chi and response partials unchanged analytically.',
                  'low_jet_effect':'Add a+b*T to(0,0); b*T to(0,1),(0,2); all u-containing low jets unchanged.'},
        'selector':{'lambda':.01,'response_importance':.1,'normalization':'Same frozen031 mean weighted-relative six-combination residual plus normalized integrated normal-fourth roughness.',
            'noninvariance':'d_tau^4(b*T)=b*T, so temperature-normal roughness penalizes the otherwise thermodynamically irrelevant entropy gauge. Potential-scale and E/TS relative weights also change; effect not attributed solely to roughness.'},
        'reference_route':'Subtract the exact known a+b*T from shifted source low jets BEFORE computing scales/fit, then restore analytically at queries. Recover comparisons in original gauge without fitting offset parameters.',
        'limits':{'max_new_fits':5,'full_grid_fits':0,'max_unknowns':1064,'no_new_fit_at_or_after_local':'19:55','threads':16,'per_fit_seconds':90},
        'fits':{}}
    products={'coarse_logQ':grid_q,'coarse_logT':grid_t,'withheld_logQ':sample_q,'withheld_logT':sample_t,
        'withheld_iq':soft_archive['withheld_iq'],'withheld_it':soft_archive['withheld_it'],
        'node_logQ':node_q,'node_logT':node_t,'archived031_coefficients':soft_archive['coarse_lambda_0.01_coefficients']}
    for field in FIELDS:
        products['printed_reference_'+field]=soft_archive['printed_reference_'+field]
        products['archived033_hard_'+field]=hard_archive['lambda_0.01_withheld_'+field]
    cases=[('baseline',0.,False),('total_b1',b,False),('total_b100',100*b,False),
           ('reference_b1',b,True),('reference_b100',100*b,True)]
    baseline=None
    for label,shift_b,remainder in cases:
        current_time=datetime.now()
        if (current_time.hour,current_time.minute)>=(19,55):
            result['fits'][label]={'status':'deferred_deadline'}
            continue
        shift_a=a if shift_b else 0.
        shifted=shifted_low(coarse_low,grid_t,shift_a,shift_b)
        fit_low=remove_input_gauge(shifted,grid_t,shift_a,shift_b) if remainder else shifted
        roundtrip={str(orders):summary((fit_low[orders]-coarse_low[orders])/
            np.maximum(1.,abs(coarse_low[orders]))) for orders in LOW_ORDERS} if remainder else None
        try:
            coefficients,solve=solve_regularized(fit_low,grid_q*LN10,grid_t*LN10,.01,response_weight=.1)
        except (ValueError,RuntimeError,TimeoutError,np.linalg.LinAlgError) as error:
            result['fits'][label]={'status':'failed','reason':str(error)}
            if label=='baseline':break
            continue
        recovered,raw,_=recovered_states(coefficients,grid_q,grid_t,sample_q,sample_t,shift_a,shift_b,remainder)
        recovered_nodes,raw_nodes,_=recovered_states(coefficients,grid_q,grid_t,node_q,node_t,shift_a,shift_b,remainder)
        if label=='baseline':baseline=recovered
        diagnostic={'status':'solved','solve':solve,'a':shift_a,'b':shift_b,'analytic_reference_before_solve':remainder,
            'removed_known_gauge_departures_from_baseline':departures(recovered,baseline),
            'removed_known_gauge_printed_holdout_errors':{field:summary((recovered[field]-soft_archive['printed_reference_'+field])/
                abs(soft_archive['printed_reference_'+field])) for field in FIELDS},
            'recovered_holdout_positivity':positive_status(recovered),'recovered_node_positivity':positive_status(recovered_nodes),
            'raw_gauge_dependent_field_status':{'heldout':positive_status(raw),'node':positive_status(raw_nodes)},
            'note_raw_signs':'Entropy/energy gauge-dependent signs are not invariant physical positivity failures; compare after known offset removal.',
            'analytic_pure_bT_roughness_with_current_fit_scale':analytic_gauge_roughness(grid_q,grid_t,shift_b,solve['potential_scale']),
            'input_subtraction_roundtrip':roundtrip}
        direct_removal={'E':raw['E']-shift_a,'S':raw['S']+shift_b}
        diagnostic['known_shift_removal_numerical_crosscheck']={field:summary((direct_removal[field]-recovered[field])/
            np.maximum(1.,abs(recovered[field]))) for field in ['E','S']}
        result['fits'][label]=diagnostic
        products[label+'_coefficients']=coefficients
        for field in OUTPUT_FIELDS:
            products[label+'_recovered_'+field]=recovered[field]
            products[label+'_raw_'+field]=raw[field]
        print(f'{label}: rank {solve["rank"]}/{solve["unknowns"]}; stationarity {solve["augmented_stationarity_relative"]:.2e}; recovered Cv p99 departure {diagnostic["removed_known_gauge_departures_from_baseline"]["Cv"]["p99"]:.3e}',flush=True)
    if 'baseline' in result['fits'] and result['fits']['baseline']['status']=='solved':
        old=soft_archive['coarse_lambda_0.01_coefficients']
        result['baseline_reproduces_frozen031_coefficients_relative_max']=float(np.max(abs(products['baseline_coefficients']-old)/np.maximum(1.,abs(old))))
    result['elapsed_seconds']=time.perf_counter()-started
    result['peak_process_rss_kib_linux']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if digest(helper)!=FROZEN_SHA:raise ValueError('Frozen helper changed during audit')
    (OUTPUT/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    np.savez_compressed(OUTPUT/'paired_gauges.npz',**products)
    lines=['# Known-gauge sensitivity of031 coarse C3 selector','',result['claim'],'',
        'Same204 active-geometry holdouts; lambda .01, response importance .1; five coarse fits only.','',
        '| Route | Removed-gauge Pgas p99 / max | Cv p99 / max | chiRho p99 / max | chiT p99 / max |',
        '| --- | --- | --- | --- | --- |']
    for label,fit in result['fits'].items():
        if fit['status']=='solved':
            stats=fit['removed_known_gauge_departures_from_baseline']
            lines.append('| '+label+' | '+' | '.join(f"{stats[field]['p99']:.3e} / {stats[field]['max']:.3e}" for field in ['Pgas','Cv','chiRho','chiT'])+' |')
        else:lines.append(f'| {label} | {fit["status"]} | - | - | - |')
    lines.extend(['','Entries are gauge-removed model departures from the baseline, not accuracy against fresh FreeEOS.',
        'Normal fourth derivative of b*T is nonzero; scale/relative E/TS weights also change. Analytic reference removal keeps original-gauge fit conventions.',
        'Raw shifted E/S signs are gauge-dependent, not invariant positivity diagnoses.'])
    (OUTPUT/'README.md').write_text('\n'.join(lines)+'\n')
    print('\n'.join(lines))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--overwrite',action='store_true')
    run(parser.parse_args())


if __name__=='__main__':main()
