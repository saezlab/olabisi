#!/usr/bin/env python

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

INDIVIDUAL_ID_MAP = {
    "SD1": "individual_01", "SD2": "individual_01", "SD3": "individual_02", "SD4": "individual_02",
    "SD5": "individual_03", "SD6": "individual_03", "SD7": "individual_04", "SD8": "individual_04",
    "SD9": "individual_05", "SD10": "individual_05", "SD11": "individual_06", "SD12": "individual_06",
    "SD13": "individual_07", "SD14": "individual_07", "SD15": "individual_08", "SD16": "individual_08",
    "SD17": "individual_09", "SD18": "individual_09", "SD19": "individual_10", "SD20": "individual_10",
    "SD21": "individual_11", "SD22": "individual_11", "SD23": "individual_12", "SD24": "individual_12",
    "SD25": "individual_13", "SD26": "individual_13", "SD27": "individual_14", "SD28": "individual_14",
    "SD29": "individual_15", "SD30": "individual_15", "SD31": "individual_16", "SD32": "individual_16",
    "SD33": "individual_17", "SD34": "individual_17", "SD35": "individual_18", "SD36": "individual_18",
    "SD37": "individual_19", "SD38": "individual_19", "SD39": "individual_20", "SD40": "individual_20",
    "SD41": "individual_21", "SD42": "individual_21", "SD43": "individual_22", "SD44": "individual_22",
    "SD45": "individual_23", "SD46": "individual_23", "SD47": "individual_24", "SD48": "individual_24",
    "SD49": "individual_25", "SD50": "individual_25", "SD51": "individual_26", "SD52": "individual_26",
    "SD53": "individual_27", "SD54": "individual_27", "SD55": "individual_28", "SD56": "individual_28",
    "SD57": "individual_29", "SD58": "individual_29", "SD59": "individual_30", "SD60": "individual_30",
    "SD61": "individual_31", "SD62": "individual_31", "SD63": "individual_32", "SD64": "individual_32",
    "SD65": "individual_33", "SD66": "individual_33",
}

root = Path(".")
in_adata = root / "output/preprocessing/objects/QC/adata.h5ad"
in_deg = root / "output/DEA/objects/pyDESeq2_interaction/deg.tsv"
out_dir = root / "output/DEA/objects/pyDESeq2_interaction_diagnostics"
out_dir.mkdir(parents=True, exist_ok=True)

adata = ad.read_h5ad(in_adata)
obs = adata.obs.copy()
obs = obs.loc[obs["QC_group"] == "A"].copy()
obs = obs.loc[
    obs["condition_group"].isin(["IFNG_neg", "IFNG_pos"])
    & obs["disease_group"].isin(["HR_control", "HR_case"])
].copy()
obs.index.name = None
obs["sample_id"] = obs.index
obs["individualID"] = obs.index.map(INDIVIDUAL_ID_MAP)
obs["pair_label"] = obs["disease_group"].astype(str) + "__" + obs["condition_group"].astype(str)

pair_counts = (
    obs.groupby(["individualID", "disease_group", "condition_group"], observed=True)
    .size()
    .rename("n_samples")
    .reset_index()
)
pair_counts.to_csv(out_dir / "pair_counts.tsv", sep="\t", index=False)

balance = (
    obs.groupby(["disease_group", "condition_group"], observed=True)
    .size()
    .rename("n_samples")
    .reset_index()
)
balance.to_csv(out_dir / "group_balance.tsv", sep="\t", index=False)

pair_completeness = (
    pair_counts.groupby("individualID", observed=True)["condition_group"]
    .nunique()
    .rename("n_conditions")
    .reset_index()
)
pair_completeness["complete_pair"] = pair_completeness["n_conditions"] == 2
pair_completeness.to_csv(out_dir / "pair_completeness.tsv", sep="\t", index=False)

deg = pd.read_csv(in_deg, sep="\t")
top_genes = deg.nsmallest(20, "pvalue")["gene_id"].tolist()

expr = adata.to_df(layer="cts.raw").loc[obs.index, top_genes].copy()
expr["sample_id"] = expr.index
expr = expr.melt(id_vars="sample_id", var_name="gene_id", value_name="count")
expr = expr.merge(
    obs[["sample_id", "individualID", "disease_group", "condition_group"]].reset_index(drop=True),
    on="sample_id",
    how="left",
)

paired = (
    expr.pivot_table(
        index=["individualID", "disease_group", "gene_id"],
        columns="condition_group",
        values="count",
        aggfunc="mean",
    )
    .reset_index()
)
if {"IFNG_neg", "IFNG_pos"}.issubset(paired.columns):
    paired["delta_ifng"] = paired["IFNG_pos"] - paired["IFNG_neg"]
else:
    paired["delta_ifng"] = np.nan
paired.to_csv(out_dir / "top20_gene_paired_deltas.tsv", sep="\t", index=False)

delta_summary = (
    paired.groupby(["disease_group", "gene_id"], observed=True)["delta_ifng"]
    .agg(["mean", "std", "median", "count"])
    .reset_index()
)
delta_summary.to_csv(out_dir / "top20_gene_delta_summary.tsv", sep="\t", index=False)
