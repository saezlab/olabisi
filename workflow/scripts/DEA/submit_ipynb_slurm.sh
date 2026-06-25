#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  ./submit_ipynb_slurm.sh [options] NOTEBOOK.ipynb [NOTEBOOK2.ipynb ...]

Submit one Slurm batch job per notebook. Each job converts the notebook to a
temporary .py script beside the notebook, runs it, and deletes the temporary
script and Slurm log when the job exits.

Defaults:
  --partition cpu-single
  --ntasks 1
  --cpus-per-task 2
  --mem 30G
  --time 2:00:00
  --conda-env squidpy_v06

Options:
  --partition NAME        Slurm partition
  --ntasks N             Number of Slurm tasks
  --cpus-per-task N      CPUs per task
  --mem SIZE             Memory request
  --time TIME            Wall time request
  --conda-env NAME       Conda environment used for nbconvert and python
  --keep-temp-py         Keep the generated temporary .py script
  --keep-log             Keep the Slurm log file
  --dry-run              Print what would be submitted
  -h, --help             Show this help

Example:
  ./submit_ipynb_slurm.sh plot_DEA.ipynb
  ./submit_ipynb_slurm.sh --cpus-per-task 4 --time 8:00:00 plot_DEA.ipynb
EOF
}

partition="cpu-single"
ntasks="1"
cpus_per_task="2"
mem="30G"
time_limit="2:00:00"
conda_env="squidpy_v06"
keep_temp_py="0"
keep_log="0"
dry_run="0"
notebooks=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --partition)
            partition="$2"
            shift 2
            ;;
        --ntasks)
            ntasks="$2"
            shift 2
            ;;
        --cpus-per-task)
            cpus_per_task="$2"
            shift 2
            ;;
        --mem)
            mem="$2"
            shift 2
            ;;
        --time)
            time_limit="$2"
            shift 2
            ;;
        --conda-env)
            conda_env="$2"
            shift 2
            ;;
        --keep-temp-py)
            keep_temp_py="1"
            shift
            ;;
        --keep-log)
            keep_log="1"
            shift
            ;;
        --dry-run)
            dry_run="1"
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        --)
            shift
            while [[ $# -gt 0 ]]; do
                notebooks+=("$1")
                shift
            done
            ;;
        -*)
            echo "Unknown option: $1" >&2
            usage >&2
            exit 2
            ;;
        *)
            notebooks+=("$1")
            shift
            ;;
    esac
done

if [[ ${#notebooks[@]} -eq 0 ]]; then
    usage >&2
    exit 2
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
project_root="$(cd "$script_dir/../../.." && pwd -P)"

for notebook in "${notebooks[@]}"; do
    if [[ "$notebook" != *.ipynb ]]; then
        echo "Expected a .ipynb file, got: $notebook" >&2
        exit 2
    fi

    if [[ ! -f "$notebook" ]]; then
        if [[ -f "$script_dir/$notebook" ]]; then
            notebook="$script_dir/$notebook"
        else
            echo "Notebook not found: $notebook" >&2
            exit 2
        fi
    fi

    notebook_dir="$(cd "$(dirname "$notebook")" && pwd -P)"
    notebook_abs="$notebook_dir/$(basename "$notebook")"
    notebook_base="$(basename "$notebook_abs" .ipynb)"
    job_name="ipynb_${notebook_base//[^A-Za-z0-9_]/_}"
    log_pattern="$notebook_dir/slurm-%j.out"

    sbatch_cmd=(
        sbatch
        --job-name "$job_name"
        --output "$log_pattern"
        --partition "$partition"
        --ntasks "$ntasks"
        --cpus-per-task "$cpus_per_task"
        --mem "$mem"
        --time "$time_limit"
        --export "ALL,PROJECT_ROOT=$project_root,NOTEBOOK=$notebook_abs,CONDA_ENV=$conda_env,KEEP_TEMP_PY=$keep_temp_py,KEEP_LOG=$keep_log"
    )

    if [[ "$dry_run" == "1" ]]; then
        printf 'Would submit:'
        printf ' %q' "${sbatch_cmd[@]}"
        printf ' <<JOB_SCRIPT\n'
        continue
    fi

    "${sbatch_cmd[@]}" <<'JOB_SCRIPT'
#!/usr/bin/env bash
set -euo pipefail

cd "$PROJECT_ROOT"

tmp_py="$(dirname "$NOTEBOOK")/.$(basename "$NOTEBOOK" .ipynb).slurm_${SLURM_JOB_ID:-$$}.py"
log_file="$(dirname "$NOTEBOOK")/slurm-${SLURM_JOB_ID:-unknown}.out"

cleanup() {
    if [[ "${KEEP_TEMP_PY:-0}" != "1" ]]; then
        rm -f "$tmp_py"
    fi
    if [[ "${KEEP_LOG:-0}" != "1" ]]; then
        rm -f "$log_file"
    fi
}
trap cleanup EXIT

echo "Started: $(date)"
echo "Node: $(hostname)"
echo "Slurm job ID: ${SLURM_JOB_ID:-unknown}"
echo "Slurm job name: ${SLURM_JOB_NAME:-unknown}"
echo "CPUs per task: ${SLURM_CPUS_PER_TASK:-unknown}"
echo "Tasks: ${SLURM_NTASKS:-unknown}"
echo "Memory per node: ${SLURM_MEM_PER_NODE:-unknown}"
echo "Working directory: $(pwd)"
echo "Project root: $PROJECT_ROOT"
echo "Notebook: $NOTEBOOK"
echo "Conda env: $CONDA_ENV"
echo "Temporary script: $tmp_py"
echo "Log file: $log_file"
echo "Keep temporary script: ${KEEP_TEMP_PY:-0}"
echo "Keep log: ${KEEP_LOG:-0}"

export PYTHONPATH="$PROJECT_ROOT${PYTHONPATH:+:$PYTHONPATH}"

echo "Python: $(conda run -n "$CONDA_ENV" python --version)"
echo "Jupyter nbconvert: $(conda run -n "$CONDA_ENV" jupyter nbconvert --version)"

conda run -n "$CONDA_ENV" jupyter nbconvert \
    --to python \
    --output-dir "$(dirname "$NOTEBOOK")" \
    --output "$(basename "$tmp_py" .py)" \
    "$NOTEBOOK"

conda run -n "$CONDA_ENV" python "$tmp_py"

echo "Finished: $(date)"
JOB_SCRIPT
done
