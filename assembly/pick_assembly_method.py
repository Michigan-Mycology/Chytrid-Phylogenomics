import numpy as np
import pandas as pd
import re
import sys

method_pat = re.compile("[.](dipspades)|[.](spades)")

def methspl(row):
    s = re.search (method_pat, row.SampleID)
    if s is not None:
        meth = s.group(1)
        if meth is None:
            meth = s.group(2)

        row["asmmeth"] = meth
        row["SampleID"]  = re.sub(method_pat, "", row["SampleID"])

    else:
        row["asmmeth"] = "NA"

    return row

def decide(x):
    if float(x["dip_hap_N50_ratio"]) > 1.0 and (float(x["dip_hap_N50_ratio"]) >= dip_hap_N50_ratio_min or float(x["dip_hap_LEN_ratio"]) <= dip_hap_LEN_ratio_max):
        x["decision"] = "dipspades"     
    else:
        x["decision"] = "spades"
    
    return x

dip_hap_N50_ratio_min = 1.10
dip_hap_LEN_ratio_max = 0.90

asmstats = pd.read_csv(sys.argv[1], sep="\t")

pert = asmstats[["SampleID", "CONTIG COUNT", "TOTAL LENGTH", "N50"]]
pert = pert[pert.SampleID.str.endswith("spades")]

pert = pert.apply(methspl, axis=1)

dips = pert[pert.asmmeth == "dipspades"].drop("asmmeth", axis=1)
spads = pert[pert.asmmeth == "spades"].drop("asmmeth", axis=1)


dips.columns = ['SampleID', 'CONTIG.COUNT_J_dip', 'TOTAL.LENGTH_J_dip', 'N50_J_dip']
spads.columns = ['SampleID', 'CONTIG.COUNT_J_spa_2', 'TOTAL.LENGTH_J_spa_2', 'N50_J_spa_2']

merge = dips.merge(spads, on="SampleID", how="left")

print(merge)
merge.to_csv("asm_stats_sub.tsv", index=False)
sys.exit()

merge = merge.assign(dip_hap_N50_ratio = lambda x: x.N50_dip/x.N50)
merge = merge.assign(dip_hap_LEN_ratio = lambda x: x["TOTAL LENGTH_dip"]/x["TOTAL LENGTH"])

merge = merge.apply(decide, axis=1)
merge = merge.assign(ID = lambda x: x.SampleID.str.split("_"))
merge["ID"] = merge["ID"].apply(lambda x: "_".join(x[2:]))

print(merge[["SampleID","decision","N50_dip", "N50", "dip_hap_N50_ratio","TOTAL LENGTH_dip","TOTAL LENGTH","dip_hap_LEN_ratio"]])

sheet = pd.read_excel("~/Downloads/2020_Chytrid_Genomes_sheets.xlsx")[["ID", "assembly_method"]]
nas_sheet = sheet[pd.isna(sheet.assembly_method)]

nas_sheet = nas_sheet.merge(merge[["ID","decision"]], on = "ID", how="left")
print(nas_sheet)


