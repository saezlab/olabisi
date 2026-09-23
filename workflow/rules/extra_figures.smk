rule plot_extra_qc:
    input:
        adata="results/preprocessing/adata.h5ad"
    output:
        pca_csv="results/extra_figures/source_data/pca.tsv",
        anova="results/extra_figures/source_data/pca_anova.tsv",
        umap_csv="results/extra_figures/source_data/umap_qc.tsv",
        pca_png="plots/extra_figures/pca_qc.png",
        pca_pdf="plots/extra_figures/pca_qc.pdf",
        umap_png="plots/extra_figures/umap_qc.png",
        umap_pdf="plots/extra_figures/umap_qc.pdf"
    params:
        pca_components=config["extra_figures"]["pca_components"]
    conda:
        "../envs/pyenv.yaml"
    script:
        "../scripts/extra_figures/plot_qc.py"


rule plot_de_diagnostics:
    input:
        degs="results/dea/interaction.tsv"
    output:
        csv="results/extra_figures/source_data/de_diagnostics.tsv",
        summary="results/extra_figures/de_significance_summary.tsv",
        pvalue_png="plots/extra_figures/de_pvalue_histogram.png",
        pvalue_pdf="plots/extra_figures/de_pvalue_histogram.pdf",
        ma_png="plots/extra_figures/de_ma.png",
        ma_pdf="plots/extra_figures/de_ma.pdf",
        volcano_png="plots/extra_figures/de_volcano.png",
        volcano_pdf="plots/extra_figures/de_volcano.pdf",
        counts_png="plots/extra_figures/de_significant_counts.png",
        counts_pdf="plots/extra_figures/de_significant_counts.pdf"
    params:
        label_top_genes=config["extra_figures"]["label_top_genes"]
    conda:
        "../envs/pyenv.yaml"
    script:
        "../scripts/extra_figures/plot_de_diagnostics.py"


rule plot_calibration_diagnostics:
    input:
        raw="results/dea/interaction.tsv",
        calibrated="results/dea/interaction_calibrated.tsv"
    output:
        csv="results/extra_figures/source_data/calibration_histograms.tsv",
        png="plots/extra_figures/calibration_histograms.png",
        pdf="plots/extra_figures/calibration_histograms.pdf"
    conda:
        "../envs/renv.yaml"
    script:
        "../scripts/extra_figures/plot_calibration_histograms.R"


rule plot_donor_diagnostics:
    input:
        adata="results/preprocessing/adata.h5ad",
        degs="results/dea/interaction_calibrated.tsv"
    output:
        deltas_csv="results/extra_figures/source_data/donor_interaction_deltas.tsv",
        report="results/extra_figures/donor_interaction_report.tsv",
        deltas_png="plots/extra_figures/donor_interaction_deltas.png",
        deltas_pdf="plots/extra_figures/donor_interaction_deltas.pdf",
        groups_csv="results/extra_figures/source_data/selected_gene_4groups.tsv",
        groups_png="plots/extra_figures/selected_gene_4groups.png",
        groups_pdf="plots/extra_figures/selected_gene_4groups.pdf",
        scatter_csv="results/extra_figures/source_data/interaction_response_scatter.tsv",
        scatter_png="plots/extra_figures/interaction_response_scatter.png",
        scatter_pdf="plots/extra_figures/interaction_response_scatter.pdf"
    params:
        selected_genes=config["figures"]["selected_genes"],
        pseudocount=config["extra_figures"]["pseudocount"],
        seed=config["random_seed"]
    conda:
        "../envs/pyenv.yaml"
    script:
        "../scripts/extra_figures/plot_donor_diagnostics.py"


rule plot_jak_stat_fgsea:
    input:
        degs="results/dea/interaction_calibrated.tsv",
        network=ancient(PROGENY_SNAPSHOT)
    output:
        csv="results/extra_figures/source_data/jak_stat_fgsea.tsv",
        png="plots/extra_figures/jak_stat_fgsea.png",
        pdf="plots/extra_figures/jak_stat_fgsea.pdf"
    params:
        seed=config["random_seed"]
    conda:
        "../envs/renv.yaml"
    script:
        "../scripts/extra_figures/plot_jak_stat_fgsea.R"
