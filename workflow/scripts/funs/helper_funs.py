from pathlib import Path


def get_project_root(caller_globals, marker="workflow/Snakefile"):
    script = caller_globals.get("__file__") or caller_globals.get("__vsc_ipynb_file__")

    if script is None:
        raise RuntimeError("Could not infer script/notebook path.")

    script = Path(script).resolve()

    for path in [script.parent, *script.parents]:
        if (path / marker).exists():
            return path

    raise RuntimeError(f"Could not find project root above {script}")
