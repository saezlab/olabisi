if not Path(ENSEMBL_SNAPSHOT).exists() or not Path(ENSEMBL_PROVENANCE).exists():
    rule retrieve_ensembl_gene_symbols:
        output:
            network=ENSEMBL_SNAPSHOT,
            summary=ENSEMBL_PROVENANCE
        params:
            resource="ensembl",
            sha256=config["prior_knowledge"]["ensembl"]["sha256"]
        conda:
            "../envs/pyenv.yaml"
        script:
            "../scripts/prior_knowledge/retrieve.py"


if not Path(COLLECTRI_SNAPSHOT).exists() or not Path(COLLECTRI_PROVENANCE).exists():
    rule retrieve_collectri:
        output:
            network=COLLECTRI_SNAPSHOT,
            summary=COLLECTRI_PROVENANCE
        params:
            resource="collectri",
            sha256=config["prior_knowledge"]["collectri"]["sha256"]
        conda:
            "../envs/pyenv.yaml"
        script:
            "../scripts/prior_knowledge/retrieve.py"


if not Path(PROGENY_SNAPSHOT).exists() or not Path(PROGENY_PROVENANCE).exists():
    rule retrieve_progeny:
        output:
            network=PROGENY_SNAPSHOT,
            summary=PROGENY_PROVENANCE
        params:
            resource="progeny",
            sha256=config["prior_knowledge"]["progeny"]["sha256"]
        conda:
            "../envs/pyenv.yaml"
        script:
            "../scripts/prior_knowledge/retrieve.py"
