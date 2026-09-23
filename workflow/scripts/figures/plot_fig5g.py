import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import configure_matplotlib, read_checked_tsv

configure_matplotlib()


if "snakemake" in locals():
    adata_path = Path(snakemake.input["adata"])
    degs_path = Path(snakemake.input["degs"])
    collectri_path = Path(snakemake.input["collectri"])
    progeny_path = Path(snakemake.input["progeny"])
    outputs = {name: Path(value) for name, value in snakemake.output.items()}
    collectri_panels = list(snakemake.params["collectri_panels"])
    progeny_panels = list(snakemake.params["progeny_panels"])
    selected_genes = list(snakemake.params["selected_genes"])
    top_targets = int(snakemake.params["top_targets"])
    seed = int(snakemake.params["seed"])
    collectri_sha256 = snakemake.params["collectri_sha256"]
    progeny_sha256 = snakemake.params["progeny_sha256"]
else:
    adata_path = Path("results/preprocessing/adata.h5ad")
    degs_path = Path("results/dea/interaction_calibrated.tsv")
    collectri_path = Path("data/prior_knowledge/collectri_20260610.tsv")
    progeny_path = Path("data/prior_knowledge/progeny_20260610.tsv")
    outputs = {}
    for key in ["IRF1", "NFKB1", "Hypoxia", "NFkB", "TNFa", "JAK_STAT", "selected_genes"]:
        outputs[f"{key}_csv"] = Path(f"results/paper/source_data/Fig5G_{key}.tsv")
        outputs[f"{key}_png"] = Path(f"plots/paper/Fig5G_{key}.png")
        outputs[f"{key}_pdf"] = Path(f"plots/paper/Fig5G_{key}.pdf")
    collectri_panels = ["IRF1", "NFKB1"]
    progeny_panels = ["Hypoxia", "NFkB", "TNFa", "JAK-STAT"]
    selected_genes = [
        "APOL1", "APOL6", "IRF1", "STAT1", "STAT2", "STAT3", "IRF2", "IRF8",
        "IRF9", "CCL2", "CXCL9", "CXCL10", "CXCL11", "ICAM1", "IDO1", "GBP1",
        "PSMB9", "TMEM140", "CTSS",
    ]
    top_targets, seed = 18, 0
    collectri_sha256 = "9fbbab30401b9e7aeb3e1c2dfbdf686e3e8aeebbd8b13cb7ff5adfb3482c3b9f"
    progeny_sha256 = "f84121b0b542ab9695221a73f1dd7b117e4c96cd39a0e0c451aeb576b436f342"

for path in outputs.values():
    path.parent.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(adata_path)
adata = adata[adata.obs["QC_group"].astype(str).eq("A")].copy()
symbols = adata.var["gene_symbol"].astype("string")
valid = symbols.notna() & symbols.str.len().gt(0)
duplicated = symbols[valid & symbols.duplicated(keep=False)].unique().tolist()
keep = (valid & ~symbols.isin(duplicated)).fillna(False).to_numpy(dtype=bool)
adata = adata[:, keep].copy()
adata.var_names = adata.var["gene_symbol"].astype(str)

expression = pd.DataFrame(
    np.asarray(adata.layers["cts.deseq2_norm"]), index=adata.obs_names, columns=adata.var_names
)
long = expression.reset_index(names="sample_id").melt(
    id_vars="sample_id", var_name="genesymbol", value_name="expression"
)
long = long.merge(
    adata.obs[["individualID", "condition_group", "disease_group"]].reset_index(),
    on="sample_id",
    how="left",
    validate="many_to_one",
)
responses = long.pivot(
    index=["individualID", "disease_group", "genesymbol"],
    columns="condition_group",
    values="expression",
).reset_index()
responses["log2_fold_change"] = np.log2(
    (responses["IFNG_pos"] + 1e-6) / (responses["IFNG_neg"] + 1e-6)
)
responses["group"] = responses["disease_group"].map(
    {
        "HR_case": "HR case, IFNG+ vs IFNG-",
        "HR_control": "HR control, IFNG+ vs IFNG-",
    }
)

degs = pd.read_csv(degs_path, sep="\t")
degs_by_symbol = degs.dropna(subset=["genesymbol"]).drop_duplicates("genesymbol").set_index("genesymbol")
collectri = read_checked_tsv(collectri_path, collectri_sha256, ["source", "target", "weight"])
progeny = read_checked_tsv(progeny_path, progeny_sha256, ["source", "target", "weight"])


def panel_key(name):
    return name.replace("-", "_")


def plot_panel(name, genes, key, title):
    genes = [gene for gene in genes if gene in expression.columns and gene in degs_by_symbol.index]
    if not genes:
        raise ValueError(f"No genes available for Figure 5G panel {name}")
    table = responses[responses["genesymbol"].isin(genes)].copy()
    table["genesymbol"] = pd.Categorical(table["genesymbol"], categories=genes, ordered=True)
    table["stat"] = table["genesymbol"].astype(str).map(degs_by_symbol["stat"])
    table["padj"] = table["genesymbol"].astype(str).map(degs_by_symbol["padj"])
    table = table.sort_values(["genesymbol", "group", "individualID"], kind="stable")
    table.to_csv(outputs[f"{key}_csv"], sep="\t", index=False)

    columns = 6
    rows = int(np.ceil(len(genes) / columns))
    fig, axes = plt.subplots(rows, columns, figsize=(2.3 * columns, 3 * rows), squeeze=False)
    palette = {
        "HR case, IFNG+ vs IFNG-": "#ff7f0e",
        "HR control, IFNG+ vs IFNG-": "#1f77b4",
    }
    rng = np.random.default_rng(seed)
    for ax, gene in zip(axes.flat, genes):
        gene_table = table[table["genesymbol"].astype(str).eq(gene)]
        order = ["HR case, IFNG+ vs IFNG-", "HR control, IFNG+ vs IFNG-"]
        sns.boxplot(
            data=gene_table,
            x="group",
            y="log2_fold_change",
            order=order,
            hue="group",
            palette=palette,
            legend=False,
            width=0.55,
            fliersize=0,
            ax=ax,
        )
        for position, group in enumerate(order):
            values = gene_table.loc[gene_table["group"].eq(group), "log2_fold_change"].to_numpy()
            jitter = rng.uniform(-0.08, 0.08, len(values))
            ax.scatter(
                np.full(len(values), position) + jitter,
                values,
                color=palette[group],
                edgecolor="black",
                linewidth=0.35,
                s=25,
                zorder=3,
            )
        stat = float(degs_by_symbol.loc[gene, "stat"])
        padj = float(degs_by_symbol.loc[gene, "padj"])
        ax.set_title(f"{gene}\nStat: {stat:.2f}\nPadj: {padj:.2e}", fontsize=8)
        ax.set_xlabel("")
        ax.set_ylabel("Log2 IFNG response")
        ax.set_xticklabels([])
    for ax in axes.flat[len(genes):]:
        ax.remove()
    fig.suptitle(title, y=1.01)
    fig.tight_layout()
    fig.savefig(outputs[f"{key}_png"], dpi=300, bbox_inches="tight")
    fig.savefig(outputs[f"{key}_pdf"], bbox_inches="tight")
    plt.close(fig)


for name in collectri_panels:
    network = collectri[collectri["source"].eq(name)].copy()
    network["stat"] = network["target"].map(degs_by_symbol["stat"])
    network["rank_score"] = network["stat"] * network["weight"]
    genes = network.sort_values("rank_score", ascending=False, kind="stable").head(top_targets)["target"].tolist()
    plot_panel(name, genes, panel_key(name), f"Top CollecTRI {name} targets")

for name in progeny_panels:
    network = progeny[progeny["source"].eq(name)].copy()
    network["stat"] = network["target"].map(degs_by_symbol["stat"])
    network["rank_score"] = network["stat"] * network["weight"]
    genes = network.sort_values("rank_score", ascending=False, kind="stable").head(top_targets)["target"].tolist()
    plot_panel(name, genes, panel_key(name), f"Top PROGENy {name} targets")

plot_panel("selected genes", selected_genes, "selected_genes", "Selected IFNG-response genes")
