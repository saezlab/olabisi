from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import configure_matplotlib

configure_matplotlib()


if "snakemake" in locals():
    adata_path = Path(snakemake.input["adata"])
    degs_path = Path(snakemake.input["degs"])
    outputs = {name: Path(value) for name, value in snakemake.output.items()}
    selected_genes = list(snakemake.params["selected_genes"])
    pseudocount = float(snakemake.params["pseudocount"])
    seed = int(snakemake.params["seed"])
else:
    adata_path = Path("results/preprocessing/adata.h5ad")
    degs_path = Path("results/dea/interaction_calibrated.tsv")
    outputs = {
        "deltas_csv": Path("results/extra_figures/source_data/donor_interaction_deltas.tsv"),
        "report": Path("results/extra_figures/donor_interaction_report.tsv"),
        "deltas_png": Path("plots/extra_figures/donor_interaction_deltas.png"),
        "deltas_pdf": Path("plots/extra_figures/donor_interaction_deltas.pdf"),
        "groups_csv": Path("results/extra_figures/source_data/selected_gene_4groups.tsv"),
        "groups_png": Path("plots/extra_figures/selected_gene_4groups.png"),
        "groups_pdf": Path("plots/extra_figures/selected_gene_4groups.pdf"),
        "scatter_csv": Path("results/extra_figures/source_data/interaction_response_scatter.tsv"),
        "scatter_png": Path("plots/extra_figures/interaction_response_scatter.png"),
        "scatter_pdf": Path("plots/extra_figures/interaction_response_scatter.pdf"),
    }
    selected_genes = [
        "APOL1", "APOL6", "IRF1", "STAT1", "STAT2", "STAT3", "IRF2", "IRF8",
        "IRF9", "CCL2", "CXCL9", "CXCL10", "CXCL11", "ICAM1", "IDO1", "GBP1",
        "PSMB9", "TMEM140", "CTSS",
    ]
    pseudocount, seed = 1e-6, 0

for path in outputs.values():
    path.parent.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(adata_path)
adata = adata[adata.obs["QC_group"].astype(str).eq("A")].copy()
symbols = adata.var["gene_symbol"].astype("string")
valid = symbols.notna() & symbols.str.len().gt(0)
duplicates = symbols[valid & symbols.duplicated(keep=False)].unique()
keep = (valid & ~symbols.isin(duplicates)).fillna(False).to_numpy(dtype=bool)
adata = adata[:, keep].copy()
adata.var_names = adata.var["gene_symbol"].astype(str)

degs = pd.read_csv(degs_path, sep="\t").dropna(subset=["genesymbol"])
degs = degs.drop_duplicates("genesymbol").set_index("genesymbol")
genes = [gene for gene in selected_genes if gene in adata.var_names and gene in degs.index]
if not genes:
    raise ValueError("None of the configured selected genes are available")

metadata = adata.obs[["individualID", "disease_group", "condition_group"]].reset_index(names="sample_id")
selected_matrix = pd.DataFrame(
    np.asarray(adata[:, genes].layers["cts.deseq2_norm_log"]),
    index=adata.obs_names,
    columns=genes,
)
groups = selected_matrix.reset_index(names="sample_id").melt(
    id_vars="sample_id", var_name="genesymbol", value_name="log_normalized_expression"
).merge(metadata, on="sample_id", validate="many_to_one")
groups["group"] = (
    groups["disease_group"].astype(str).str.replace("HR_", "", regex=False)
    + " / "
    + groups["condition_group"].astype(str)
)
groups["genesymbol"] = pd.Categorical(groups["genesymbol"], categories=genes, ordered=True)
groups = groups.sort_values(["genesymbol", "group", "individualID"], kind="stable")
groups.to_csv(outputs["groups_csv"], sep="\t", index=False)

order = ["control / IFNG_neg", "control / IFNG_pos", "case / IFNG_neg", "case / IFNG_pos"]
columns = 5
rows = int(np.ceil(len(genes) / columns))
fig, axes = plt.subplots(rows, columns, figsize=(3.0 * columns, 2.8 * rows), squeeze=False)
for ax, gene in zip(axes.flat, genes):
    table = groups[groups["genesymbol"].astype(str).eq(gene)]
    sns.boxplot(data=table, x="group", y="log_normalized_expression", order=order, fliersize=0, ax=ax)
    sns.stripplot(
        data=table, x="group", y="log_normalized_expression", order=order,
        color="black", alpha=0.55, size=2.5, jitter=0.18, ax=ax,
    )
    ax.set_title(gene)
    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=55, labelsize=8)
for ax in axes.flat[len(genes):]:
    ax.remove()
fig.tight_layout()
fig.savefig(outputs["groups_png"], dpi=300, bbox_inches="tight")
fig.savefig(outputs["groups_pdf"], bbox_inches="tight")
plt.close(fig)

expression = pd.DataFrame(
    np.asarray(adata.layers["cts.deseq2_norm"]), index=adata.obs_names, columns=adata.var_names
)
long = expression.reset_index(names="sample_id").melt(
    id_vars="sample_id", var_name="genesymbol", value_name="expression"
).merge(metadata, on="sample_id", validate="many_to_one")
responses = long.pivot(
    index=["individualID", "disease_group", "genesymbol"],
    columns="condition_group",
    values="expression",
).reset_index()
responses["log2_ifng_response"] = np.log2(
    (responses["IFNG_pos"] + pseudocount) / (responses["IFNG_neg"] + pseudocount)
)

deltas = responses[responses["genesymbol"].isin(genes)].copy()
deltas["group_mean_response"] = deltas.groupby(
    ["genesymbol", "disease_group"], observed=True
)["log2_ifng_response"].transform("mean")
deltas["centered_response"] = deltas["log2_ifng_response"] - deltas["group_mean_response"]
deltas["interaction_stat"] = deltas["genesymbol"].map(degs["stat"])
deltas["interaction_padj"] = deltas["genesymbol"].map(degs["padj"])
deltas = deltas.sort_values(["disease_group", "individualID", "genesymbol"], kind="stable")
deltas.to_csv(outputs["deltas_csv"], sep="\t", index=False)

report = deltas.pivot_table(
    index="genesymbol",
    columns="disease_group",
    values="log2_ifng_response",
    aggfunc=["mean", "std", "count"],
    observed=True,
)
report.columns = [f"{stat}_{group}" for stat, group in report.columns]
report = report.reset_index()
report["mean_interaction_delta"] = report["mean_HR_case"] - report["mean_HR_control"]
report["deseq2_interaction_stat"] = report["genesymbol"].map(degs["stat"])
report["deseq2_interaction_padj"] = report["genesymbol"].map(degs["padj"])
report.sort_values("deseq2_interaction_stat", ascending=False).to_csv(
    outputs["report"], sep="\t", index=False
)

heatmap = deltas.pivot(index=["disease_group", "individualID"], columns="genesymbol", values="log2_ifng_response")
heatmap = heatmap.reindex(columns=genes)
fig, ax = plt.subplots(figsize=(max(9, len(genes) * 0.55), 9))
sns.heatmap(heatmap, center=0, cmap="vlag", ax=ax, cbar_kws={"label": "Log2 IFNG response"})
ax.set(xlabel="Gene", ylabel="Disease group / donor", title="Per-donor IFNG-response diagnostics")
fig.tight_layout()
fig.savefig(outputs["deltas_png"], dpi=300, bbox_inches="tight")
fig.savefig(outputs["deltas_pdf"], bbox_inches="tight")
plt.close(fig)

scatter = responses.groupby(["genesymbol", "disease_group"], observed=True)["log2_ifng_response"].mean().unstack()
scatter = scatter.rename(columns={"HR_case": "case_mean_response", "HR_control": "control_mean_response"}).reset_index()
scatter["interaction_delta"] = scatter["case_mean_response"] - scatter["control_mean_response"]
scatter["interaction_stat"] = scatter["genesymbol"].map(degs["stat"])
scatter["interaction_padj"] = scatter["genesymbol"].map(degs["padj"])
scatter["selected_gene"] = scatter["genesymbol"].isin(genes)
scatter = scatter.sort_values(["selected_gene", "interaction_stat", "genesymbol"], ascending=[False, False, True])
scatter.to_csv(outputs["scatter_csv"], sep="\t", index=False)

fig, ax = plt.subplots(figsize=(7, 7))
rng = np.random.default_rng(seed)
ax.scatter(scatter["control_mean_response"], scatter["case_mean_response"], s=8, alpha=0.25, color="#808080")
selected = scatter[scatter["selected_gene"]]
ax.scatter(selected["control_mean_response"], selected["case_mean_response"], s=30, color="#d62728")
for row in selected.itertuples():
    ax.annotate(row.genesymbol, (row.control_mean_response, row.case_mean_response), fontsize=8,
                xytext=(rng.uniform(2, 6), rng.uniform(2, 6)), textcoords="offset points")
limits = [
    min(scatter["control_mean_response"].min(), scatter["case_mean_response"].min()),
    max(scatter["control_mean_response"].max(), scatter["case_mean_response"].max()),
]
ax.plot(limits, limits, linestyle="--", color="black", linewidth=0.8)
ax.set(
    xlabel="Mean log2 IFNG response in controls",
    ylabel="Mean log2 IFNG response in cases",
    title="Case-control interaction-response comparison",
)
fig.tight_layout()
fig.savefig(outputs["scatter_png"], dpi=300, bbox_inches="tight")
fig.savefig(outputs["scatter_pdf"], bbox_inches="tight")
plt.close(fig)
