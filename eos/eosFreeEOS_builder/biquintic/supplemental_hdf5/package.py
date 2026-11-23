"""Package only supplements and small metadata, with deterministic tar headers."""
import argparse
import json
from pathlib import Path
import subprocess
import tarfile
import sys
sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from freeeos_potential import checksum
from supplemental_hdf5.verify import verify


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir",type=Path,required=True)
    parser.add_argument("--supplement-dir",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--audit-output",type=Path,required=True)
    parser.add_argument("--threads",type=int,default=8,choices=range(1,9))
    parser.add_argument("--extraction-dir",type=Path,required=True)
    args=parser.parse_args()
    result=verify(args.source_dir,args.supplement_dir)
    members=sorted(list(args.supplement_dir.glob("*.bq.h5"))+list(args.supplement_dir.glob("*.bq.audit.json"))+
                   [args.supplement_dir/"supplements.manifest",args.supplement_dir/"supplements.audit.json"],key=lambda p:p.name)
    assert len(members)==180 and all(p.is_file() and not p.is_symlink() for p in members)
    tarpath=args.output.with_name(args.output.name+".tar.partial")
    partial=args.output.with_name(args.output.name+".partial")
    with tarfile.open(tarpath,"w",format=tarfile.GNU_FORMAT) as tar:
        for path in members:
            info=tarfile.TarInfo("eosFreeEOS_data/"+path.name)
            info.size=path.stat().st_size; info.mode=0o644; info.mtime=0
            info.uid=info.gid=0; info.uname=info.gname=""
            with path.open("rb") as stream: tar.addfile(info,stream)
    with partial.open("wb") as stream:
        subprocess.run(["xz",f"-T{args.threads}","-6","--stdout",str(tarpath)],stdout=stream,check=True)
    partial.chmod(0o644); partial.replace(args.output)
    subprocess.run(["xz","--test",str(args.output)],check=True)
    args.extraction_dir.mkdir(parents=True,exist_ok=False)
    original_dir=args.extraction_dir/"eosFreeEOS_data"
    original_dir.mkdir()
    for source in args.source_dir.glob("*.data"):
        (original_dir/source.name).symlink_to(source.resolve())
    with tarfile.open(args.output,"r:xz") as tar:
        entries=tar.getmembers()
        assert [m.name for m in entries]==["eosFreeEOS_data/"+p.name for p in members]
        assert all(m.isfile() and m.mode==0o644 and m.uid==m.gid==m.mtime==0 for m in entries)
        assert all(len(Path(m.name).parts)==2 and not m.name.endswith(".data") for m in entries)
        tar.extractall(args.extraction_dir)
    assert len(list(original_dir.glob("*.data")))==89
    extracted=verify(args.source_dir,original_dir)
    result.update({"archive_path":str(args.output.resolve()),"archive_bytes":args.output.stat().st_size,
                   "archive_MiB":args.output.stat().st_size/2**20,"archive_sha256":checksum(args.output),
                   "archive_members":len(members),"archive_hdf5_count":89,"archive_source_data_count":0,
                   "tar_bytes":tarpath.stat().st_size,"all_member_file_bytes":sum(p.stat().st_size for p in members),
                   "staged_extraction":str(args.extraction_dir.resolve()),"extraction_verified":extracted==verify(args.source_dir,args.supplement_dir)})
    tarpath.unlink()
    args.audit_output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
if __name__ == "__main__": main()
