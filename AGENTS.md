## Snakemake

- Use the standard script interface patterns:
  - R scripts should use `if (exists("snakemake")) { ... } else { ... }`.
  - Python scripts should use `if "snakemake" in locals(): ... else: ...`.
  - In the `if` branch, read all `snakemake.input`, `snakemake.output`, and `snakemake.params`.
  - In the `else` branch, keep local defaults for direct script execution.
- For Python scripts run via Snakemake `script:`, do **not** use `from __future__ import ...`.
  - Snakemake prepends runtime code, so future imports may no longer be the first statement and can fail with `SyntaxError`.
- When suggesting workflow execution, always provide exact copy-paste Snakemake commands.
- Prefer explicit targets and include key flags: `--cores` and `--software-deployment-method apptainer` or `--software-deployment-method conda`.
- Command format defaults (switch out `apptainer` and `conda` according to project preference):
  - Dry-run: `snakemake --cores 1 <target> --software-deployment-method apptainer -n`
  - Run: `snakemake --cores 1 <target> --software-deployment-method apptainer`
- Do not use the cache/home wrapper variant by default in suggested commands.
- Only provide the cache/home wrapper variant if the user explicitly asks for it or if plain `snakemake ...` fails in this environment.
- Keep Snakemake naming simple when reuse is likely:
  - Prefer generic I/O key names like `rds`, `csv`, `plot`, `degs_plot`, `top_markers_plot`, `manifest_csv`.
  - Avoid branch-specific key names unless the rule/script truly needs multiple variants of the same type.
  - Reusable scripts should accept simple keys first (for example `rds`) and only keep legacy aliases temporarily for backward compatibility.
- Keep Snakemake rule files easy to read:
  - Prefer direct literal paths in `input:` and `output:` when a path is owned by a single rule.
  - Do not create top-of-file path constants for one-off outputs just to avoid repeating a literal once.
  - Centralize only config-backed values, shared wildcard expansions, or aggregates reused across multiple rules or across a rule file and `Snakefile`.
  - Avoid: `FOO_PLOT = "plots/geomx/foo.pdf"` when only one rule uses it.
  - Prefer: `output: plot="plots/geomx/foo.pdf"` directly in the owning rule.
- For dataset-specific workflows, prefer canonical literal manuscript-facing paths over local path-root aliases.
  - Prefer `results/lupus/processing/...` or `plots/lupus/processing/...` directly in the owning rule.
  - Avoid introducing aliases like `LUPUS_RESULTS_ROOT`, `LUPUS_PLOTS_ROOT`, `ANNOTATION_ROOT`, or similar when they mainly hide the actual path layout.
  - Only centralize a path-like value when it is reused in a way that clearly improves readability.
- Keep `params:` limited to true behavioral controls.
  - Pass values like `resolution`, `resolution_grid`, `dims_use`, `group_col`, `cluster_col`, or `label_col` through `params:` when the script needs them as behavior.
  - Do not pass path roots or output directories through `params:` when the script can infer them from `snakemake.output`.
  - Avoid pass-through params such as `results_root`, `plots_root`, `qc_dir`, `embeddings_dir`, `clustering_dir`, or `markers_dir` unless they are genuinely required inputs.
- Treat `config.yaml` as real configuration, not as a mirror of fixed workflow paths.
  - Do not store dataset-specific path aliases in config when there is already one canonical manuscript-facing location.
  - Remove config keys once they no longer affect workflow behavior.
  - Prefer compact readable config forms when the script parser supports them cleanly, for example `dims_use: "1:30"` instead of a long explicit integer list.
- In Snakemake script fallback `else` blocks, prefer direct literal paths for fixed manuscript-facing defaults.
  - Prefer `annotation_rds <- "results/lupus/processing/annotation/objects/lupus_annotation_subclustered.rds"` over verbose `file.path(...)` construction when the path is static.
  - Use `file.path(...)` when paths are genuinely being composed from changing parts, not just to spell a fixed string.
- Remove stale default concepts when they no longer drive behavior.
  - Do not keep settings like `default_resolution` in config, rules, or scripts unless they change outputs, dependencies, or analysis behavior.
  - If all resolutions are processed uniformly, derive behavior from `resolution_grid` alone.
- For reusable rules/scripts, keep interfaces minimal and explicit:
  - Inputs: only the files required by that script.
  - Params: only behavioral controls (for example `cluster_column_prefix`, `resolution`, `markers_group`).
  - Outputs: stable generic names that map cleanly to file paths.

Python Snakemake skeleton:
```python
if "snakemake" in locals():
    in_path = snakemake.input["in_path"]
    out_path = snakemake.output["out_path"]
else:
    in_path = "path/to/local/input"
    out_path = "path/to/local/output"
```

R Snakemake skeleton:
```r
if (exists("snakemake")) {
  in_path <- snakemake@input[["in_path"]]
  out_path <- snakemake@output[["out_path"]]
} else {
  in_path <- "path/to/local/input"
  out_path <- "path/to/local/output"
}
```

Reusable rule/script naming example:
```python
rule example_dotplots:
    input:
        rds="results/lupus/object.rds",
        markers_csv="results/lupus/findallmarkers_res_0.6.csv",
        degs_csv="data/reference/e1_degs.csv"
    output:
        degs_plot="plots/lupus/dotplot_degs_res_0.6.pdf",
        top_markers_plot="plots/lupus/dotplot_top3_marker_genes_res_0.6.pdf",
        manifest_csv="results/lupus/dotplot_marker_manifest_res_0.6.csv"
    params:
        cluster_column_prefix="MLC_clusters_res.",
        resolution="0.6",
        markers_group="MLC"
```