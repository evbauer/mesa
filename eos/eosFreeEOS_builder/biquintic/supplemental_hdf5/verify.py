"""Validate all 89 source identities, HDF5 schemas, grids and numeric policies."""
import argparse
import json
from pathlib import Path
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import h5py
import numpy as np
from freeeos_potential import checksum
from supplemental_hdf5.produce import MIXED, POLICY_SHA, digest_bytes


def verify(source_dir, supplement_dir):
    lines = (supplement_dir / "supplements.manifest").read_text().splitlines()
    assert lines[0] == "FREEEOS_BIQUINTIC_SUPPLEMENTS_V1"
    assert int(lines[1]) == len(lines)-2 == 89
    sources = {path.stem for path in source_dir.glob("mesa-FreeEOS_*.data")}
    assert len(sources) == 89
    covered = set()
    byte_count = raw_count = co_count = nodes = cells = 0
    for line in lines[2:]:
        x, z, available, name, source_sha, payload_sha = line.split()
        assert available == "1" and name.endswith(".bq.h5")
        stem = name[:-6]
        assert stem in sources and stem not in covered
        covered.add(stem)
        source = source_dir / (stem+".data")
        supplement = supplement_dir / name
        assert checksum(source) == source_sha
        assert checksum(supplement) == payload_sha
        audit = json.loads((supplement_dir / (stem+".bq.audit.json")).read_text())
        with h5py.File(supplement, "r") as h5:
            assert set(h5.keys()) == set(MIXED) | {"valid", "source_sha256", "preprocessing_sha256"}
            assert set(h5.attrs) == set(audit["metadata"])
            for key, value in audit["metadata"].items():
                assert h5.attrs[key] == value, key
                assert h5.attrs.get_id(key).shape == ()
                assert h5.attrs.get_id(key).dtype == (np.dtype("int32") if isinstance(value,int) else np.dtype("float64"))
            assert h5.attrs["X"] == float(x) and h5.attrs["Z"] == float(z)
            assert h5.attrs["family_id"] == int(stem.endswith("_CO"))
            co_count += int(stem.endswith("_CO"))
            shape = (int(h5.attrs["ntau"]),int(h5.attrs["nu"]))
            valid = h5["valid"][:]
            assert valid.shape == shape and h5["valid"].dtype == np.dtype("int8")
            assert np.isin(valid,[0,1]).all() and valid.any()
            assert int(valid.sum()) == audit["complete_jet_nodes"]
            nodes += int(valid.sum()); cells += audit["complete_cells"]
            for key in MIXED:
                dataset = h5[key]
                assert dataset.shape == shape and dataset.dtype == np.dtype("float64")
                assert dataset.scaleoffset is None
                values = dataset[:]
                assert np.isfinite(values).all() and np.all(values[valid == 0] == 0)
            for key, sha in (("source_sha256",source_sha),("preprocessing_sha256",POLICY_SHA)):
                assert h5[key].shape == (32,) and h5[key].dtype == np.dtype("int32")
                assert np.array_equal(h5[key][:],digest_bytes(sha))
            raw_count += sum(h5[key].size*h5[key].dtype.itemsize for key in h5)
        byte_count += supplement.stat().st_size
    assert covered == sources and co_count == 6
    assert set(path.name for path in supplement_dir.glob("*.h5")) == {stem+".bq.h5" for stem in sources}
    return {"source_count":89,"supplement_count":89,"co_count":co_count,
            "supplement_file_bytes":byte_count,"supplement_file_MiB":byte_count/2**20,
            "unfiltered_dataset_bytes":raw_count,"unfiltered_dataset_MiB":raw_count/2**20,
            "complete_nodes":nodes,"complete_cells":cells,"external_sha_verification":True}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir",type=Path,required=True)
    parser.add_argument("--supplement-dir",type=Path,required=True)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    result=verify(args.source_dir,args.supplement_dir)
    content=json.dumps(result,indent=2)+"\n"
    if args.output: args.output.write_text(content)
    print(content)
if __name__ == "__main__": main()
