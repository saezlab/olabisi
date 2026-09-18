import platform
import subprocess
import sys
from pathlib import Path


if "snakemake" in locals():
    output_path = Path(snakemake.output["summary"])
else:
    output_path = Path("results/provenance/pyenv_packages.txt")

output_path.parent.mkdir(parents=True, exist_ok=True)
freeze = subprocess.run(
    [sys.executable, "-m", "pip", "freeze"], check=True, capture_output=True, text=True
).stdout
output_path.write_text(
    f"python={platform.python_version()}\nexecutable={sys.executable}\n" + freeze
)
