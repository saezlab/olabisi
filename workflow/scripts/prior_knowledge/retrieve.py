import hashlib
import json
import platform
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import decoupler as dc


if "snakemake" in locals():
    snapshot_path = Path(snakemake.output["network"])
    provenance_path = Path(snakemake.output["summary"])
    resource = snakemake.params["resource"]
    expected_sha256 = snakemake.params["sha256"]
else:
    snapshot_path = Path("data/prior_knowledge/collectri_20260610.tsv")
    provenance_path = Path("data/prior_knowledge/collectri_20260610.provenance.json")
    resource = "collectri"
    expected_sha256 = "9fbbab30401b9e7aeb3e1c2dfbdf686e3e8aeebbd8b13cb7ff5adfb3482c3b9f"

snapshot_path.parent.mkdir(parents=True, exist_ok=True)
provenance_path.parent.mkdir(parents=True, exist_ok=True)

if resource == "ensembl":
    arguments = {
        "dataset": "hsapiens_gene_ensembl",
        "host": "http://www.ensembl.org",
        "attributes": [
            "ensembl_gene_id",
            "external_gene_name",
            "gene_biotype",
            "chromosome_name",
        ],
    }
    required = ["gene_id", "gene_symbol", "gene_biotype", "chromosome"]
    provider = "Ensembl BioMart"
elif resource == "collectri":
    arguments = {"organism": "human", "remove_complexes": False, "license": "academic"}
    required = ["source", "target", "weight"]
    provider = "decoupler/OmniPath"
elif resource == "progeny":
    arguments = {"organism": "human", "top": 500, "thr_padj": 0.05}
    required = ["source", "target", "weight", "padj"]
    provider = "decoupler/OmniPath"
else:
    raise ValueError(f"Unknown prior-knowledge resource: {resource}")

if snapshot_path.exists():
    import pandas as pd

    table = pd.read_csv(snapshot_path, sep="\t")
    acquisition = "validated_existing_snapshot"
elif resource == "ensembl":
    from pybiomart import Dataset

    dataset = Dataset(name=arguments["dataset"], host=arguments["host"])
    table = dataset.query(attributes=arguments["attributes"]).rename(
        columns={
            "Gene stable ID": "gene_id",
            "Gene name": "gene_symbol",
            "Gene type": "gene_biotype",
            "Chromosome/scaffold name": "chromosome",
        }
    )
    table = table.drop_duplicates("gene_id")
    acquisition = "retrieved"
elif resource == "collectri":
    table = dc.op.collectri(**arguments)[["source", "target", "weight"]]
    acquisition = "retrieved"
else:
    table = dc.op.progeny(**arguments)
    acquisition = "retrieved"

missing = set(required).difference(table.columns)
if missing:
    raise ValueError(f"Retrieved {resource} table is missing columns: {sorted(missing)}")
if acquisition == "retrieved":
    table.to_csv(snapshot_path, sep="\t", index=False)
digest = hashlib.sha256(snapshot_path.read_bytes()).hexdigest()
if digest != expected_sha256:
    raise ValueError(
        f"Retrieved {resource} checksum {digest} does not match frozen snapshot {expected_sha256}. "
        "Restore the timestamped Zenodo snapshot instead of silently changing the analysis."
    )

provenance = {
    "resource": resource,
    "acquisition": acquisition,
    "retrieved_at": datetime.now(timezone.utc).isoformat(),
    "provider": provider,
    "arguments": arguments,
    "schema": required,
    "sha256": digest,
    "python": platform.python_version(),
    "packages": {"decoupler": version("decoupler")},
}
provenance_path.write_text(json.dumps(provenance, indent=2) + "\n")
