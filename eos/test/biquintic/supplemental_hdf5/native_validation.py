"""Compare native source/HDF5 against explicit V1, suffixes, handles and toggles."""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys
sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from runtime_bridge.runtime_environment import private_environment,installed_inventory


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--supplement-dir",type=Path,required=True)
    parser.add_argument("--reference-dir",type=Path,required=True)
    args=parser.parse_args()
    output=args.output.resolve()
    workspace=Path(__file__).resolve().parents[5]
    mesa=workspace/"mesa"
    massfile=mesa/"eos/eosFreeEOS_builder/mass_frac_solar.txt"
    before=installed_inventory()
    metrics={}
    queries=[]
    for rho,t in [(-7,5),(-5,5.3),(-3,5.7),(-1,6.1)]:
        for hydrogen,metals in [(.7,.02),(.65,.03),(.3,.1),(0,.75),(0,.8),(.05,.85),(0,.95),(0,1)]:
            queries.append(("component",rho,t,hydrogen,metals))

    def run_case(name,directory=args.supplement_dir,enabled=True,suffix="configured",threads=1,old=False,query_list=queries):
        case=output/"native_api"/name
        case.mkdir(parents=True,exist_ok=True)
        environment=private_environment(case)
        environment["OMP_NUM_THREADS"]=str(threads)
        controls=["&eos"," use_WDM_AV=.false.",f" FreeEOS_use_HELM_interpolation={'.true.' if enabled else '.false.'}"]
        if directory is not None: controls.append(f" FreeEOS_biquintic_table_dir='{directory.resolve()}'")
        if suffix=="ordinary": controls.extend(f" suffix_for_FreeEOS_Z({layer})=''" for layer in (13,14,15))
        if suffix=="mixed": controls.append(" suffix_for_FreeEOS_Z(13)=''")
        inlist=case/"inlist_eos"
        inlist.write_text("\n".join(controls+["/"])+"\n")
        inputs="\n".join(" ".join(map(str,row)) for row in query_list)+"\n"
        (case/"queries.txt").write_text(inputs)
        driver=output/"tests/native_driver"
        if old: driver=workspace/"freeeos_biquintic_progress/experiments/075_co_table_support/native_build/driver/native_driver"
        command=[str(driver),str(mesa),str(inlist),str(massfile)]
        (case/"command.json").write_text(json.dumps(command,indent=2)+"\n")
        run=subprocess.run(command,input=inputs,capture_output=True,text=True,env=environment,timeout=300)
        (case/"stdout.txt").write_text(run.stdout); (case/"stderr.txt").write_text(run.stderr)
        (case/"returncode.txt").write_text(str(run.returncode)+"\n")
        if run.returncode: raise RuntimeError(name+": "+run.stderr[-2000:])
        rows=np.asarray([[float(value) for value in line.split()] for line in run.stdout.splitlines() if line[:1].isdigit()])
        assert len(rows)==len(query_list)
        return rows

    for suffix in ["configured","ordinary","mixed"]:
        current=run_case("hdf5_"+suffix,suffix=suffix)
        reference=run_case("v1_"+suffix,directory=args.reference_dir,suffix=suffix)
        assert np.array_equal(current[:,1],reference[:,1])
        successful=current[:,1]==0
        assert successful.all(), (suffix,current[:,1].tolist())
        actual=current[successful,5:]; expected=reference[successful,5:]
        assert np.array_equal(np.isnan(actual),np.isnan(expected))
        residual=np.abs(actual-expected)/np.maximum(1,np.abs(expected))
        metrics[suffix]={"query_count":len(current),"statuses":current[:,1].astype(int).tolist(),
                         "max_scaled_output_difference":float(np.nanmax(residual)),
                         "max_scaled_per_output":np.nanmax(residual,axis=0).tolist(),
                         "table_loads":int(current[-1,4])}
        if suffix=="configured": configured=current
        if suffix=="ordinary": ordinary=current
    threaded=run_case("hdf5_threads4",threads=4)
    assert np.array_equal(threaded,configured,equal_nan=True)
    metrics["native_thread_count_exact_match"]=True
    default=run_case("installed_default",directory=None)
    assert np.array_equal(default,configured,equal_nan=True)
    metrics["installed_default_exact_match"]=True
    off=run_case("toggle_off",enabled=False)
    previous=run_case("previous_toggle_off",enabled=False,old=True)
    assert np.array_equal(off[:,1],previous[:,1])
    assert np.array_equal(off[:,5:],previous[:,5:],equal_nan=True)
    metrics["toggle_off_previous_exact_match"]=True
    case=output/"native_api/dual_handles"
    case.mkdir(parents=True,exist_ok=True)
    command=[str(output/"tests/dual_handle_driver"),str(mesa),
             str(output/"native_api/hdf5_configured/inlist_eos"),str(output/"native_api/hdf5_ordinary/inlist_eos"),str(massfile)]
    run=subprocess.run(command,capture_output=True,text=True,env=private_environment(case),timeout=300)
    (case/"command.json").write_text(json.dumps(command,indent=2)+"\n")
    (case/"stdout.txt").write_text(run.stdout); (case/"stderr.txt").write_text(run.stderr)
    assert run.returncode==0
    rows=[line.split() for line in run.stdout.splitlines() if line[:1].isdigit()]
    assert len(rows)==6
    # Query index 5 is X=.05/Z=.85; the dual driver uses X=0/Z=.85, so obtain its reference explicitly.
    selected_dual=run_case("dual_reference_configured",query_list=[("component",-7,5,0,.85)])
    ordinary_dual=run_case("dual_reference_ordinary",suffix="ordinary",query_list=[("component",-7,5,0,.85)])
    for index,row in enumerate(rows):
        assert int(row[1])==0
        expected=selected_dual[0,5] if index%2==0 else ordinary_dual[0,5]
        assert float(row[2])==expected
    assert "COUNTS 6 0 4" in run.stdout
    metrics["dual_handle_counts"]=[6,0,4]
    assert installed_inventory()==before
    metrics["installed_library_and_DT_cache_unchanged"]=True
    (output/"native_api/results.json").write_text(json.dumps(metrics,indent=2)+"\n")
    print(json.dumps({key:value for key,value in metrics.items() if not isinstance(value,dict)},indent=2))
if __name__=="__main__": main()
