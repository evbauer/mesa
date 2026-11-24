"""Python/FORUM asymmetric fixture and recoverable schema/identity errors."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import h5py
import numpy as np
from freeeos_potential import CRAD
from supplemental_hdf5.produce import MIXED,POLICY_SHA,digest_bytes


def write_fixture(path):
    attrs={"format_version":1,"source_version":51,"family_id":0,"nu":7,"ntau":5,
           "radiation_split":1,"stencil_radius":2,"twist_route":1,"preprocessing_version":1,
           "X":.7,"Z":.02,"u_min":-2.0,"du":.1,"tau_min":8.0,"dtau":.2,
           "radiation_constant":float(CRAD),"minimum_log10_T":3.2}
    iu,it=np.meshgrid(np.arange(1,8),np.arange(1,6))
    valid=(iu+2*it)%4!=0
    with h5py.File(path,"w") as h5:
        for key,value in attrs.items():
            h5.attrs[key]=np.int32(value) if isinstance(value,int) else np.float64(value)
        for index,key in enumerate(MIXED,1):
            h5[key]=np.where(valid,100000*index+10*iu+it+.125,0).astype(np.float64)
        h5["valid"]=valid.astype(np.int8)
        h5["source_sha256"]=digest_bytes("a"*64)
        h5["preprocessing_sha256"]=digest_bytes(POLICY_SHA)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--driver",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--sdk",type=Path,required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    environment=os.environ.copy()
    environment["LD_LIBRARY_PATH"]=f"{args.sdk.resolve()}/lib:{args.sdk.resolve()}/lib64:"+environment.get("LD_LIBRARY_PATH","")
    environment.update(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
    cases={"good":0,"missing_attribute":11,"rank_attribute":11,"wrong_family":12,"wrong_composition":12,
           "wrong_version":11,"wrong_radius":11,"wrong_floor":11,"wrong_reference":11,"wrong_digest":12,
           "digest_range":11,"digest_rank":11,"wrong_policy":12,"missing_dataset":11,"nonsquare_transpose":11,
           "rank_dataset":11,"invalid_mask":11,"nan_accepted":11,"nonzero_invalid":11,"oversize":11,"missing_file":10}
    results=[]
    for name,status in cases.items():
        path=args.output/(name+".h5")
        if name!="missing_file":
            write_fixture(path)
            with h5py.File(path,"r+") as h5:
                if name=="missing_attribute": del h5.attrs["nu"]
                if name=="rank_attribute": h5.attrs["nu"]=np.array([7],dtype=np.int32)
                for case,key,value in [("wrong_family","family_id",1),("wrong_composition","X",.6),
                     ("wrong_version","format_version",2),("wrong_radius","stencil_radius",3),
                     ("wrong_floor","minimum_log10_T",3.1),("wrong_reference","radiation_constant",CRAD*1.01),
                     ("oversize","nu",100001)]:
                    if name==case: h5.attrs[key]=value
                if name=="wrong_digest": h5["source_sha256"][0]=3
                if name=="digest_range": h5["source_sha256"][0]=256
                if name=="wrong_policy": h5["preprocessing_sha256"][0]=0
                if name=="digest_rank":
                    del h5["source_sha256"]; h5["source_sha256"]=np.zeros((4,8),dtype=np.int32)
                if name=="missing_dataset": del h5["f_uu_tau"]
                if name=="nonsquare_transpose":
                    values=h5["f_uu_tau"][:].T
                    del h5["f_uu_tau"]; h5["f_uu_tau"]=values
                if name=="rank_dataset":
                    del h5["f_uu_tau"]; h5["f_uu_tau"]=np.zeros(35)
                if name=="invalid_mask": h5["valid"][0,0]=2
                if name=="nan_accepted": h5["f_uu_tau"][0,0]=np.nan
                if name=="nonzero_invalid": h5["f_uu_tau"][0,1]=1
        command=[str(args.driver.resolve()),str(path.resolve()),str(status)]
        run=subprocess.run(command,capture_output=True,text=True,env=environment,timeout=30)
        (args.output/(name+".log")).write_text(json.dumps(command)+"\n"+run.stdout+run.stderr)
        if run.returncode: raise RuntimeError(name+": "+run.stdout+run.stderr)
        results.append({"name":name,"expected_status":status,"returncode":run.returncode})
    (args.output/"results.json").write_text(json.dumps(results,indent=2)+"\n")
    print(f"PASS {len(results)} Python/FORUM fixture cases")
if __name__=="__main__": main()
