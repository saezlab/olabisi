import hashlib
from pathlib import Path

import pandas as pd
import pytest
import yaml


ROOT = Path(__file__).resolve().parents[1]
CONFIG = yaml.safe_load((ROOT / "config/config.yaml").read_text())


def test_raw_metadata_pairing():
    metadata = pd.read_csv(ROOT / "data/metadata.csv")
    assert len(metadata) == 33
    assert metadata[["IFNG_neg", "IFNG_pos"]].notna().all().all()
    assert metadata[["IFNG_neg", "IFNG_pos"]].stack().nunique() == 66
    assert (metadata["QC_group"] == "Group A").sum() == 25


def test_prior_knowledge_snapshots_match_config():
    expected_columns = {
        "ensembl": {"gene_id", "gene_symbol", "gene_biotype", "chromosome"},
        "collectri": {"source", "target", "weight"},
        "progeny": {"source", "target", "weight", "padj"},
    }
    for name, expected in expected_columns.items():
        entry = CONFIG["prior_knowledge"][name]
        path = ROOT / entry["path"]
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        assert digest == entry["sha256"]
        assert expected.issubset(pd.read_csv(path, sep="\t", nrows=1).columns)


def test_fixed_paper_panels():
    figures = CONFIG["figures"]
    assert figures["collectri_panels"] == ["IRF1", "NFKB1"]
    assert figures["progeny_panels"] == ["Hypoxia", "NFkB", "TNFa", "JAK-STAT"]
    assert "PI3K" not in figures["progeny_panels"]
    assert figures["top_targets"] == 18


def test_python_snakemake_script_contract():
    scripts = [
        ROOT / "workflow/scripts/preprocessing/clean_metadata.py",
        ROOT / "workflow/scripts/preprocessing/make_adata.py",
        ROOT / "workflow/scripts/preprocessing/preprocess.py",
        ROOT / "workflow/scripts/DEA/pyDESeq2_interaction.py",
        ROOT / "workflow/scripts/DEA/run_ulm.py",
        ROOT / "workflow/scripts/figures/plot_fig5bc.py",
        ROOT / "workflow/scripts/figures/plot_enrichment.py",
        ROOT / "workflow/scripts/figures/plot_fig5g.py",
    ]
    for script in scripts:
        text = script.read_text()
        assert 'if "snakemake" in locals():' in text
        assert "from __future__" not in text


def test_editable_pdf_typography_contract():
    style = (ROOT / "workflow/scripts/common.py").read_text()
    assert '"font.family": "Arial"' in style
    assert '"font.size": 8' in style
    assert '"pdf.fonttype": 42' in style


@pytest.mark.skipif(
    not (ROOT / "results/enrichment/ulm.tsv").exists(),
    reason="Generated ULM table is validated after the paper workflow runs",
)
def test_generated_ulm_schema_and_panels():
    table = pd.read_csv(ROOT / "results/enrichment/ulm.tsv", sep="\t")
    assert list(table.columns) == [
        "comparison",
        "source",
        "activity",
        "padj",
        "gene_set",
        "significant_padj_0.05",
    ]
    assert "pvalue" not in table.columns
    significant = table[table["significant_padj_0.05"]]
    assert {"IRF1", "NFKB1"}.issubset(
        set(significant.loc[significant["gene_set"] == "collectri", "source"])
    )
    assert {"Hypoxia", "NFkB", "TNFa", "JAK-STAT", "PI3K"}.issubset(
        set(significant.loc[significant["gene_set"] == "progeny", "source"])
    )


@pytest.mark.skipif(
    not (ROOT / "results/dea/interaction.tsv").exists(),
    reason="Generated DE tables are validated after the paper workflow runs",
)
def test_generated_interaction_results():
    raw = pd.read_csv(ROOT / "results/dea/interaction.tsv", sep="\t")
    calibrated = pd.read_csv(
        ROOT / "results/dea/interaction_calibrated.tsv", sep="\t"
    )
    assert len(raw) == 19_556
    assert len(calibrated) == 19_556
    assert raw["padj"].lt(0.05).sum() == 0
    assert calibrated["padj"].lt(0.05).sum() == 0


@pytest.mark.skipif(
    not (ROOT / "results/preprocessing/qc_summary.tsv").exists(),
    reason="Generated QC summary is validated after the paper workflow runs",
)
def test_generated_sample_counts():
    summary = pd.read_csv(
        ROOT / "results/preprocessing/qc_summary.tsv", sep="\t"
    ).set_index("metric")["value"]
    assert int(summary["samples_total"]) == 66
    assert int(summary["samples_qc_group_a"]) == 50
    assert int(summary["complete_qc_group_a_pairs"]) == 25
