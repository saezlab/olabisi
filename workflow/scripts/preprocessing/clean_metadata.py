from pathlib import Path

import pandas as pd


if "snakemake" in locals():
    metadata_path = Path(snakemake.input["metadata"])
    output_path = Path(snakemake.output["csv"])
else:
    metadata_path = Path("data/metadata.csv")
    output_path = Path("results/preprocessing/metadata.tsv")

output_path.parent.mkdir(parents=True, exist_ok=True)
metadata = pd.read_csv(metadata_path)
required = {"disease_group", "IFNG_neg", "IFNG_pos", "QC_group"}
missing = required.difference(metadata.columns)
if missing:
    raise ValueError(f"Metadata is missing columns: {sorted(missing)}")

# A donor is one original wide row. Assign it before melting so the pairing is
# explicit and does not rely on sample-name ordering.
metadata = metadata.copy()
metadata["individualID"] = [f"individual_{index:02d}" for index in range(1, len(metadata) + 1)]
metadata = metadata.melt(
    id_vars=["individualID", "disease_group", "QC_group"],
    value_vars=["IFNG_neg", "IFNG_pos"],
    var_name="condition_group",
    value_name="sample_id",
)
metadata["disease_group"] = metadata["disease_group"].replace(
    {"HR Case": "HR_case", "HR Control": "HR_control"}
)
metadata["QC_group"] = metadata["QC_group"].replace({"Group A": "A", "Group B": "B"})
if len(metadata) != 66:
    raise ValueError(f"Expected 66 samples, found {len(metadata)}")
if metadata["sample_id"].duplicated().any():
    raise ValueError("Sample identifiers are not unique")
pair_sizes = metadata.groupby("individualID")["condition_group"].nunique()
if not pair_sizes.eq(2).all():
    raise ValueError("Every donor must have both IFNG conditions in the raw metadata")

metadata.to_csv(output_path, sep="\t", index=False)
