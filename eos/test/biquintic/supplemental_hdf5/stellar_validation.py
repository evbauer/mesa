"""Private stock ns_c bounded toggle/V1 comparison and declared photo restart."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from freeeos_potential import checksum
from runtime_bridge.runtime_environment import private_environment,installed_inventory


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--supplement-dir",type=Path,required=True)
    parser.add_argument("--reference-dir",type=Path,required=True)
    parser.add_argument("--models",type=int,default=40)
    parser.add_argument("--timeout",type=int,default=900)
    args=parser.parse_args()
    workspace=Path(__file__).resolve().parents[5]
    mesa=workspace/"mesa"; sdk=workspace/"mesasdk-x86_64-linux-26.6.1"
    output=args.output.resolve(); stellar=output/"stellar"
    stellar.mkdir(parents=True,exist_ok=True)
    build=stellar/"driver"
    build.mkdir(exist_ok=True)
    stock=mesa/"star/test_suite/ns_c"
    environment=private_environment(stellar)
    environment["PATH"]=str(sdk/"bin")+":"+environment["PATH"]
    environment["PKG_CONFIG_PATH"]=":".join(str(path) for path in (output/"native_build").glob("*/lib/pkgconfig"))+":"+str(sdk/"lib/pkgconfig")
    includes=subprocess.check_output(["pkg-config","--define-variable=math_slot=crmath","--cflags","mesa-star"],text=True,env=environment).split()
    libraries=subprocess.check_output(["pkg-config","--define-variable=math_slot=crmath","--libs","--static","mesa-star"],text=True,env=environment).split()
    command=[str(sdk/"bin/gfortran"),"-O2","-g","-fopenmp","-ftrampoline-impl=heap","-ffp-contract=off","-ffree-line-length-none","-fcheck=all",
             "-I"+str(build),"-J"+str(build),"-I"+str(mesa/"star/job"),"-I"+str(mesa/"star/private"),"-I"+str(mesa/"include"),*includes,
             str(stock/"src/run_star_extras.f90"),str(mesa/"star/job/run_star.f90"),str(stock/"src/run.f90"),"-o",str(build/"star"),
             "-L"+str(sdk/"lib"),"-L"+str(sdk/"lib64"),"-Wl,-rpath,"+str(sdk/"lib"),"-Wl,-rpath,"+str(sdk/"lib64"),*libraries]
    (build/"command.json").write_text(json.dumps(command,indent=2)+"\n")
    with (build/"build.log").open("w") as log:
        subprocess.run(command,env=environment,stdout=log,stderr=subprocess.STDOUT,check=True)
    before=installed_inventory()
    results={"executable_sha256":checksum(build/"star"),"stock_source":str(stock),"model_limit":args.models,"threads":1,"cases":{}}

    def prepare(name,enabled,directory):
        run_directory=stellar/name
        run_directory.mkdir(exist_ok=False)
        for filename in ["ns_c_env.mod","history_columns.list","profile_columns.list","inlist_to_c_flash","inlist_pgstar"]:
            shutil.copyfile(stock/filename,run_directory/filename)
        # Stock physical controls are retained; the cap and output intervals define this validation.
        source=(run_directory/"inlist_to_c_flash").read_text()
        source=source.replace("required_termination_code_string = 'power_nuc_burn_upper_limit'","required_termination_code_string = 'max_model_number'")
        (run_directory/"inlist_to_c_flash").write_text(source)
        lines=[]
        for namelist in ["star_job","eos","kap","controls","pgstar"]:
            lines.extend(["&"+namelist,f" read_extra_{namelist}_inlist(1)=.true.",f" extra_{namelist}_inlist_name(1)='inlist_to_c_flash'",
                          f" read_extra_{namelist}_inlist(2)=.true.",f" extra_{namelist}_inlist_name(2)='validation_controls'","/"])
        (run_directory/"inlist").write_text("\n".join(lines)+"\n")
        controls=["&star_job"," pgstar_flag=.false.","/","&eos"," use_WDM_AV=.false.",
                  f" FreeEOS_use_HELM_interpolation={'.true.' if enabled else '.false.'}"]
        if directory is not None: controls.append(f" FreeEOS_biquintic_table_dir='{directory.resolve()}'")
        controls.extend(["/","&kap","/","&controls",f" max_model_number={args.models}"," photo_interval=20"," profile_interval=20",
                         " history_interval=1"," terminal_interval=10","/","&pgstar","/"])
        (run_directory/"validation_controls").write_text("\n".join(controls)+"\n")
        return run_directory

    def run(name,run_directory,restart=None):
        run_environment=private_environment(run_directory)
        run_environment["OMP_NUM_THREADS"]="1"
        if restart is not None: shutil.copyfile(restart,run_directory/"restart_photo")
        command=[str(build/"star")]
        (run_directory/"command.json").write_text(json.dumps(command)+"\n")
        with (run_directory/"run.log").open("w") as log:
            process=subprocess.run(command,cwd=run_directory,env=run_environment,stdout=log,stderr=subprocess.STDOUT,timeout=args.timeout)
        text=(run_directory/"run.log").read_text()
        if process.returncode or "termination code: max_model_number" not in text:
            raise RuntimeError(f"{name} did not reach declared model cap; inspect {run_directory}/run.log")
        history=np.genfromtxt(run_directory/"LOGS/history.data",names=True,skip_header=5)
        history=np.atleast_1d(history)
        photos=sorted((run_directory/"photos").glob("x*"))
        result={"returncode":process.returncode,"last_model":int(history["model_number"][-1]),"history_rows":len(history),
                "last_age":float(history["star_age"][-1]),"photo_names":[path.name for path in photos],
                "final_model_sha256":checksum(run_directory/"final.mod"),"termination":"max_model_number"}
        results["cases"][name]=result
        (stellar/"checkpoint.json").write_text(json.dumps(results,indent=2)+"\n")
        print(name,json.dumps(result),flush=True)
        return history

    histories={}
    for name,enabled,directory in [("hdf5",True,None),("v1",True,args.reference_dir),("toggle_off",False,None)]:
        run_directory=prepare(name,enabled,directory)
        histories[name]=run(name,run_directory)
    for name in ["hdf5","toggle_off"]:
        restart_directory=prepare(name+"_restart",name=="hdf5",None)
        checkpoint=stellar/name/"photos/x020"
        if not checkpoint.exists(): raise RuntimeError("Declared x020 photo missing")
        histories[name+"_restart"]=run(name+"_restart",restart_directory,checkpoint.resolve())
        original=histories[name]; resumed=histories[name+"_restart"]
        metrics={}
        for column in original.dtype.names:
            if resumed[column].dtype.kind not in "fiu": continue
            matches=np.isin(original["model_number"],resumed["model_number"])
            left=original[column][matches]
            right=resumed[column][np.isin(resumed["model_number"],original["model_number"][matches])]
            if len(left)!=len(right): raise RuntimeError("Restart model alignment")
            metrics[column]=float(np.nanmax(np.abs(left-right)/np.maximum(1,np.abs(left))))
        results[name+"_restart_max_scaled_history_difference"]=max(metrics.values())
        results[name+"_restart_per_column"]=metrics
    metrics={}
    for column in histories["hdf5"].dtype.names:
        if histories["hdf5"][column].dtype.kind in "fiu":
            left=histories["hdf5"][column]; right=histories["v1"][column]
            assert len(left)==len(right)
            metrics[column]=float(np.nanmax(np.abs(left-right)/np.maximum(1,np.abs(right))))
    results["hdf5_v1_per_column_scaled_history_difference"]=metrics
    assert installed_inventory()==before
    results["installed_library_and_DT_cache_unchanged"]=True
    (stellar/"results.json").write_text(json.dumps(results,indent=2)+"\n")
    print("PASS bounded stellar toggle/V1 and x020 restart validation",flush=True)
if __name__=="__main__": main()
