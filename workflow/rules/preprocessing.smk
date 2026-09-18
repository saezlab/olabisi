rule clean_metadata:
    input:
        metadata="data/metadata.csv"
    output:
        csv="results/preprocessing/metadata.tsv"
    conda:
        "../envs/pyenv.yaml"
    script:
        "../scripts/preprocessing/clean_metadata.py"


rule make_adata:
    input:
        metadata="results/preprocessing/metadata.tsv",
        counts="data/gene_counts.txt",
        network=ancient(ENSEMBL_SNAPSHOT)
    output:
        adata="results/preprocessing/adata_unfiltered.h5ad"
    params:
        sha256=config["prior_knowledge"]["ensembl"]["sha256"]
    conda:
        "../envs/pyenv.yaml"
    script:
        "../scripts/preprocessing/make_adata.py"


rule preprocess:
    input:
        adata="results/preprocessing/adata_unfiltered.h5ad"
    output:
        adata="results/preprocessing/adata.h5ad",
        counts="results/preprocessing/counts.tsv",
        summary="results/preprocessing/qc_summary.tsv"
    params:
        seed=config["random_seed"],
        min_cells=config["preprocessing"]["min_cells"],
        min_total_count=config["preprocessing"]["min_total_count"],
        min_count=config["preprocessing"]["min_count"],
        min_prop=config["preprocessing"]["min_prop"],
        highly_variable_genes=config["preprocessing"]["highly_variable_genes"],
        pca_components=config["preprocessing"]["pca_components"],
        neighbors=config["preprocessing"]["neighbors"],
        leiden_resolution=config["preprocessing"]["leiden_resolution"]
    resources:
        mem_mb=32000
    conda:
        "../envs/pyenv.yaml"
    script:
        "../scripts/preprocessing/preprocess.py"
