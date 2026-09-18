rule record_pyenv_packages:
    output:
        summary="results/provenance/pyenv_packages.txt"
    conda:
        "../envs/pyenv.yaml"
    script:
        "../scripts/provenance/python_manifest.py"


rule record_renv_packages:
    output:
        summary="results/provenance/renv_packages.txt"
    conda:
        "../envs/renv.yaml"
    script:
        "../scripts/provenance/r_manifest.R"


rule comparison_report:
    input:
        old_degs="archive/20260918T143651Z_pre_refactor/output/DEA/objects/pyDESeq2_interaction/deg.tsv",
        new_degs="results/dea/interaction.tsv",
        old_calibrated="archive/20260918T143651Z_pre_refactor/output/DEA/objects/calibrate_pvalues_interaction/deg_calibrated.tsv",
        new_calibrated="results/dea/interaction_calibrated.tsv",
        old_fig5b="archive/20260918T143651Z_pre_refactor/output/preprocessing/plots/QC/heatmap.png",
        new_fig5b="plots/paper/Fig5B_heatmap.png",
        old_fig5c_groups="archive/20260918T143651Z_pre_refactor/output/preprocessing/plots/QC/umap_02.png",
        new_fig5c_groups="plots/paper/Fig5C_umap_groups.png",
        old_fig5c_markers="archive/20260918T143651Z_pre_refactor/output/preprocessing/plots/QC/umap_04.png",
        new_fig5c_markers="plots/paper/Fig5C_umap_markers.png"
    output:
        summary="results/provenance/old_to_new_comparison.json"
    params:
        rtol=1e-6,
        atol=1e-8
    conda:
        "../envs/pyenv.yaml"
    script:
        "../scripts/provenance/compare_archive.py"
