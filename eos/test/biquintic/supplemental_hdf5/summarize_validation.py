"""Summarize all-table error scales, populations and response-derivative outliers."""
import argparse
import json
from pathlib import Path
import sys
sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from freeeos_potential import CRAD,LN10
from family_cache.produce import read_audited

STATE_NAMES=["P","Pgas","E","S","Cv","Cp","chiRho","chiT","gamma1","gamma3","grad_ad"]
NODE_NAMES=["f","f_u","f_uu","f_tau","f_u_tau","f_uu_tau","f_tau_tau","f_u_tau_tau","f_uu_tau_tau"]
PARTIAL_NAMES=["f","f_x","f_t","f_xx","f_xt","f_tt","f_xxx","f_xxt","f_xtt","f_ttt"]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log",type=Path,required=True)
    parser.add_argument("--source-dir",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    tables=[]; partials=[]; nodes=[]
    for line in args.log.read_text().splitlines():
        columns=line.split()
        if not columns: continue
        if columns[0]=="TABLE":
            values=np.asarray(columns[5:],dtype=float)
            assert len(values)==82
            tables.append({"filename":columns[1],"samples":int(columns[2]),"accepted_states":int(columns[3]),"nonphysical_states":int(columns[4]),"values":values})
        if columns[0]=="PARTIAL_OUTLIER":
            context=list(map(float,columns[4:]))
            partials.append({"filename":columns[1],"field":STATE_NAMES[int(columns[2])-1],"axis":"lnrho" if columns[3]=="1" else "lnT",
                             **dict(zip(["log10_rho","log10_T","actual","reference","absolute_difference","scaled_difference","relative_difference"],context))})
        if columns[0]=="NODE_OUTLIER":
            context=list(map(float,columns[3:]))
            nodes.append({"filename":columns[1],"constraint":NODE_NAMES[int(columns[2])-1],
                          **dict(zip(["iq_one_based","it_one_based","actual","reference","absolute_difference","scaled_difference"],context))})
    assert len(tables)==len(partials)==len(nodes)==89
    assert "PASS all source+HDF5/V1 masks and mixed constraints: 89" in args.log.read_text()
    maxima=np.max([record["values"] for record in tables],axis=0)
    groups={}; start=0
    for name,names in [("nodal_constraints",NODE_NAMES),("potential_partials",PARTIAL_NAMES),("state",STATE_NAMES),("state_partials",STATE_NAMES)]:
        groups[name]={"names":names,"maximum_scaled_difference":maxima[start:start+len(names)].tolist(),
                      "maximum_absolute_difference":maxima[start+41:start+41+len(names)].tolist()}
        start+=len(names)
    worst=max(partials,key=lambda record:record["scaled_difference"])
    source=args.source_dir/(worst["filename"][:-6]+".data")
    header,logq,logt,rows,audit=read_audited(source)
    q=worst["log10_rho"]-2*worst["log10_T"]+12
    iq=int(np.searchsorted(logq,q)-1); it=int(np.searchsorted(logt,worst["log10_T"])-1)
    rho=10**worst["log10_rho"]; temperature=10**worst["log10_T"]
    radiation_pressure=CRAD*temperature**4/3
    corners=[]
    for i in [iq,iq+1]:
        for j in [it,it+1]:
            row=rows[i,j]; corner_temperature=10**logt[j]; corner_density=10**(logq[i]+2*logt[j]-12)
            gas_pressure=10**row[1]; prad=CRAD*corner_temperature**4/3
            energy=10**row[2]; entropy=10**row[3]
            corners.append({"iq_zero_based":i,"it_zero_based":j,"log10_Q":float(logq[i]),"log10_T":float(logt[j]),
                            "Pgas":float(gas_pressure),"Prad":float(prad),"Prad_over_Pgas":float(prad/gas_pressure),
                            "E":float(energy),"T_times_S":float(corner_temperature*entropy),"radiation_scale":float(prad/corner_density),
                            "Cv":float(row[7]),"chiRho":float(row[4]),"chiT":float(row[5]),
                            "f_after_split":float(energy-corner_temperature*entropy+prad/corner_density)})
    result={"table_count":89,"sampled_cells":sum(record["samples"] for record in tables),
            "accepted_states":sum(record["accepted_states"] for record in tables),"nonphysical_states":sum(record["nonphysical_states"] for record in tables),
            "population":"All accepted nodal jets; about64 deterministic interior-cell queries per table at offsets(.61,.37). This is not exhaustive cell sampling.",
            "scaled_difference_definition":"abs(actual-reference)/max(1,abs(reference)), computed componentwise. Absolute maxima and scaled maxima need not occur at the same point.",
            "masks_and_mixed_exact":True,"groups":groups,"worst_state_partial":worst,
            "worst_state_partial_source_context":{"source":source.name,"source_sha256":audit["source_sha256"],"log10_Q":q,
                                                 "rho":rho,"T":temperature,"Prad":radiation_pressure,"corners":corners},
            "worst_nodal_scaled_difference":max(nodes,key=lambda record:record["scaled_difference"]),
            "partial_outliers_by_table":partials,"nodal_outliers_by_table":nodes}
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({key:value for key,value in result.items() if key not in ["partial_outliers_by_table","nodal_outliers_by_table","groups"]},indent=2))
if __name__=="__main__": main()
