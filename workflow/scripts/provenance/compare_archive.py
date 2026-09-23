import json
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
from skimage.metrics import structural_similarity


if "snakemake" in locals():
    paths = {name: Path(value) for name, value in snakemake.input.items()}
    output_path = Path(snakemake.output["summary"])
    rtol = float(snakemake.params["rtol"])
    atol = float(snakemake.params["atol"])
else:
    paths = {
        "old_degs": Path("archive/20260918T143651Z_pre_refactor/output/DEA/objects/pyDESeq2_interaction/deg.tsv"),
        "new_degs": Path("results/dea/interaction.tsv"),
        "old_calibrated": Path("archive/20260918T143651Z_pre_refactor/output/DEA/objects/calibrate_pvalues_interaction/deg_calibrated.tsv"),
        "new_calibrated": Path("results/dea/interaction_calibrated.tsv"),
        "old_fig5b": Path("archive/20260918T143651Z_pre_refactor/output/preprocessing/plots/QC/heatmap.png"),
        "new_fig5b": Path("plots/paper/Fig5B_heatmap.png"),
        "old_fig5c_groups": Path("archive/20260918T143651Z_pre_refactor/output/preprocessing/plots/QC/umap_02.png"),
        "new_fig5c_groups": Path("plots/paper/Fig5C_umap_groups.png"),
        "old_fig5c_markers": Path("archive/20260918T143651Z_pre_refactor/output/preprocessing/plots/QC/umap_04.png"),
        "new_fig5c_markers": Path("plots/paper/Fig5C_umap_markers.png"),
    }
    output_path = Path("results/provenance/old_to_new_comparison.json")
    rtol, atol = 1e-6, 1e-8

output_path.parent.mkdir(parents=True, exist_ok=True)


def compare_table(old_path, new_path, key, numeric_columns):
    old = pd.read_csv(old_path, sep="\t")
    new = pd.read_csv(new_path, sep="\t")
    if key not in old.columns:
        old = old.rename(columns={old.columns[0]: key})
    merged = old[[key, *[column for column in numeric_columns if column in old.columns]]].merge(
        new[[key, *[column for column in numeric_columns if column in new.columns]]],
        on=key,
        how="outer",
        suffixes=("_old", "_new"),
        indicator=True,
    )
    report = {
        "old_rows": int(len(old)),
        "new_rows": int(len(new)),
        "matched_rows": int(merged["_merge"].eq("both").sum()),
        "only_old": int(merged["_merge"].eq("left_only").sum()),
        "only_new": int(merged["_merge"].eq("right_only").sum()),
        "columns": {},
    }
    for column in numeric_columns:
        old_column, new_column = f"{column}_old", f"{column}_new"
        if old_column not in merged or new_column not in merged:
            report["columns"][column] = {"status": "not_comparable"}
            continue
        both = merged[merged["_merge"].eq("both")][[old_column, new_column]].dropna()
        differences = np.abs(both[old_column].to_numpy() - both[new_column].to_numpy())
        report["columns"][column] = {
            "compared": int(len(both)),
            "allclose": bool(np.allclose(both[old_column], both[new_column], rtol=rtol, atol=atol)),
            "max_abs_difference": float(differences.max()) if len(differences) else None,
        }
    return report


def compare_image(old_path, new_path):
    old = np.asarray(Image.open(old_path).convert("RGB"))
    new = np.asarray(Image.open(new_path).convert("RGB"))
    report = {"old_dimensions": list(old.shape), "new_dimensions": list(new.shape)}
    if old.shape == new.shape:
        report["structural_similarity"] = float(
            structural_similarity(old, new, channel_axis=2, data_range=255)
        )
    else:
        report["structural_similarity"] = None
        report["note"] = "Not computed because dimensions differ"
    return report


report = {
    "tolerances": {"rtol": rtol, "atol": atol},
    "tables": {
        "interaction": compare_table(
            paths["old_degs"],
            paths["new_degs"],
            "gene_id",
            ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"],
        ),
        "interaction_calibrated": compare_table(
            paths["old_calibrated"],
            paths["new_calibrated"],
            "gene_id",
            ["stat", "pvalue", "padj", "pvalue_raw", "padj_raw", "null_sd"],
        ),
        "ulm": {
            "status": "no_archived_reference",
            "note": "origin/main ULM outputs were not present in the pre-refactor output tree",
        },
    },
    "images": {
        "Fig5B": compare_image(paths["old_fig5b"], paths["new_fig5b"]),
        "Fig5C_groups": compare_image(paths["old_fig5c_groups"], paths["new_fig5c_groups"]),
        "Fig5C_markers": compare_image(paths["old_fig5c_markers"], paths["new_fig5c_markers"]),
    },
    "authority": "TSV agreement is authoritative; image comparison is structural only",
}
output_path.write_text(json.dumps(report, indent=2) + "\n")
