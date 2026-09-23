rule paired_interaction_deseq2:
    input:
        adata="results/preprocessing/adata.h5ad"
    output:
        degs="results/dea/interaction.tsv"
    params:
        seed=config["random_seed"]
    resources:
        mem_mb=32000
    conda:
        "../envs/pyenv.yaml"
    script:
        "../scripts/DEA/pyDESeq2_interaction.py"


rule calibrate_interaction_pvalues:
    input:
        degs="results/dea/interaction.tsv"
    output:
        degs="results/dea/interaction_calibrated.tsv",
        summary="results/dea/calibration_summary.tsv"
    conda:
        "../envs/renv.yaml"
    script:
        "../scripts/DEA/calibrate_pvalues_interaction.R"


rule run_ulm:
    input:
        degs="results/dea/interaction_calibrated.tsv",
        collectri=ancient(COLLECTRI_SNAPSHOT),
        progeny=ancient(PROGENY_SNAPSHOT)
    output:
        csv="results/enrichment/ulm.tsv",
        degs="results/enrichment/ulm_input.tsv",
        summary="results/enrichment/ulm_summary.tsv"
    params:
        collectri_sha256=config["prior_knowledge"]["collectri"]["sha256"],
        progeny_sha256=config["prior_knowledge"]["progeny"]["sha256"]
    conda:
        "../envs/pyenv.yaml"
    script:
        "../scripts/DEA/run_ulm.py"
