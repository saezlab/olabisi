import sys
from pathlib import Path

import decoupler as dc
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import read_checked_tsv


if "snakemake" in locals():
    degs_path = Path(snakemake.input["degs"])
    collectri_path = Path(snakemake.input["collectri"])
    progeny_path = Path(snakemake.input["progeny"])
    output_path = Path(snakemake.output["csv"])
    input_table_path = Path(snakemake.output["degs"])
    summary_path = Path(snakemake.output["summary"])
    collectri_sha256 = snakemake.params["collectri_sha256"]
    progeny_sha256 = snakemake.params["progeny_sha256"]
else:
    degs_path = Path("results/dea/interaction_calibrated.tsv")
    collectri_path = Path("data/prior_knowledge/collectri_20260610.tsv")
    progeny_path = Path("data/prior_knowledge/progeny_20260610.tsv")
    output_path = Path("results/enrichment/ulm.tsv")
    input_table_path = Path("results/enrichment/ulm_input.tsv")
    summary_path = Path("results/enrichment/ulm_summary.tsv")
    collectri_sha256 = "9fbbab30401b9e7aeb3e1c2dfbdf686e3e8aeebbd8b13cb7ff5adfb3482c3b9f"
    progeny_sha256 = "f84121b0b542ab9695221a73f1dd7b117e4c96cd39a0e0c451aeb576b436f342"

for path in [output_path, input_table_path, summary_path]:
    path.parent.mkdir(parents=True, exist_ok=True)

degs = pd.read_csv(degs_path, sep="\t")
required = {"gene_id", "genesymbol", "stat", "contrast"}
missing = required.difference(degs.columns)
if missing:
    raise ValueError(f"DE table is missing columns: {sorted(missing)}")
rows_input = len(degs)
degs = degs.dropna(subset=["genesymbol", "stat"]).copy()
degs = degs[degs["genesymbol"].astype(str).str.len().gt(0)].copy()
duplicated_symbols = sorted(
    degs.loc[degs["genesymbol"].duplicated(keep=False), "genesymbol"].unique().tolist()
)
degs = degs[~degs["genesymbol"].isin(duplicated_symbols)].copy()
degs = degs.sort_values("gene_id", kind="stable")
degs.to_csv(input_table_path, sep="\t", index=False)

matrix = degs.pivot(index="contrast", columns="genesymbol", values="stat")
collectri = read_checked_tsv(collectri_path, collectri_sha256, ["source", "target", "weight"])
progeny = read_checked_tsv(progeny_path, progeny_sha256, ["source", "target", "weight"])

results = []
for gene_set, network in [("collectri", collectri), ("progeny", progeny)]:
    activities, adjusted_pvalues = dc.mt.ulm(data=matrix, net=network, verbose=False)
    activity_long = activities.rename_axis("comparison").reset_index().melt(
        id_vars="comparison", var_name="source", value_name="activity"
    )
    padj_long = adjusted_pvalues.rename_axis("comparison").reset_index().melt(
        id_vars="comparison", var_name="source", value_name="padj"
    )
    table = activity_long.merge(padj_long, on=["comparison", "source"], validate="one_to_one")
    table["gene_set"] = gene_set
    table["significant_padj_0.05"] = table["padj"].lt(0.05)
    results.append(table)

enrichment = pd.concat(results, ignore_index=True)
enrichment = enrichment.sort_values(["gene_set", "padj", "source"], kind="stable").reset_index(drop=True)
enrichment.to_csv(output_path, sep="\t", index=False)

summary = pd.DataFrame(
    {
        "metric": [
            "de_rows_input", "de_rows_ulm", "duplicated_symbols_removed",
            "collectri_sources", "progeny_sources",
            "collectri_bh_significant", "progeny_bh_significant",
        ],
        "value": [
            rows_input, len(degs), len(duplicated_symbols),
            int(enrichment["gene_set"].eq("collectri").sum()),
            int(enrichment["gene_set"].eq("progeny").sum()),
            int((enrichment["gene_set"].eq("collectri") & enrichment["significant_padj_0.05"]).sum()),
            int((enrichment["gene_set"].eq("progeny") & enrichment["significant_padj_0.05"]).sum()),
        ],
    }
)
summary.to_csv(summary_path, sep="\t", index=False)
