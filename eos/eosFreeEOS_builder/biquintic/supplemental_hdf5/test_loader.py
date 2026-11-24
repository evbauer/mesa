"""Historical raw-reader authentication/grid tests; not production loader checks."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import h5py
from freeeos_potential import checksum
from supplemental_hdf5.produce import digest_bytes


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--driver",type=Path,required=True)
    parser.add_argument("--source-dir",type=Path,required=True)
    parser.add_argument("--supplement-dir",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--sdk",type=Path,required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    tokens=(args.supplement_dir/"supplements.manifest").read_text().splitlines()[2].split()
    hydrogen,metals,available,filename,source_sha,supplement_sha=tokens
    original_source=args.source_dir/(filename[:-6]+".data")
    original_supplement=args.supplement_dir/filename
    cases={"good":0,"wrong_supplement_hash":12,"changed_source":12,"changed_grid_metadata":12,"changed_mask":12,
           "authenticated_wrong_source_composition":12,"missing_source":10,"missing_supplement":10}
    environment=dict(os.environ)
    environment["LD_LIBRARY_PATH"]=f"{args.sdk.resolve()}/lib:{args.sdk.resolve()}/lib64:"+environment.get("LD_LIBRARY_PATH","")
    environment.update(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
    records=[]
    for name,status in cases.items():
        source=original_source; supplement=original_supplement; expected_source=source_sha; expected_supplement=supplement_sha
        if name=="wrong_supplement_hash": expected_supplement="1"*64
        if name=="changed_source":
            source=args.output/"changed_source.data"
            shutil.copyfile(original_source,source)
            with source.open("ab") as stream: stream.write(b"\n")
        if name in ["changed_grid_metadata","changed_mask","authenticated_wrong_source_composition"]:
            supplement=args.output/(name+".h5")
            shutil.copyfile(original_supplement,supplement)
            with h5py.File(supplement,"r+") as h5:
                if name=="changed_grid_metadata": h5.attrs["u_min"]+=.01
                if name=="changed_mask": h5["valid"][0,0]=1
                if name=="authenticated_wrong_source_composition":
                    source=args.output/"wrong_composition.data"
                    with original_source.open("rb") as src,source.open("wb") as dst:
                        dst.write(src.readline())
                        header=src.readline().split(); header[1]=b"0.1000"
                        dst.write(b" ".join(header)+b"\n")
                        shutil.copyfileobj(src,dst)
                    expected_source=checksum(source)
                    h5["source_sha256"][:]=digest_bytes(expected_source)
            expected_supplement=checksum(supplement)
        if name=="missing_source": source=args.output/"missing.data"
        if name=="missing_supplement": supplement=args.output/"missing.h5"
        command=[str(args.driver.resolve()),str(source.resolve()),str(supplement.resolve()),hydrogen,metals,"0",expected_source,expected_supplement,str(status)]
        run=subprocess.run(command,capture_output=True,text=True,env=environment,timeout=60)
        (args.output/(name+".log")).write_text(json.dumps(command)+"\n"+run.stdout+run.stderr)
        if run.returncode: raise RuntimeError(name+": "+run.stdout+run.stderr)
        records.append({"name":name,"expected_status":status,"returncode":run.returncode})
    (args.output/"results.json").write_text(json.dumps(records,indent=2)+"\n")
    print(f"PASS {len(records)} source+HDF5 loader/authentication cases")
if __name__=="__main__": main()
