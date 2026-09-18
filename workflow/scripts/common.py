import hashlib
from pathlib import Path

import pandas as pd


def ensure_parent(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_checked_tsv(path, expected_sha256, required_columns):
    observed = sha256(path)
    if observed != expected_sha256:
        raise ValueError(
            f"Checksum mismatch for {path}: expected {expected_sha256}, observed {observed}"
        )
    table = pd.read_csv(path, sep="\t")
    missing = set(required_columns).difference(table.columns)
    if missing:
        raise ValueError(f"{path} is missing required columns: {sorted(missing)}")
    return table
