from pathlib import Path


def get_project_root(caller_globals, marker="workflow/Snakefile"):
    location = (
        caller_globals.get("__file__")
        or caller_globals.get("__vsc_ipynb_file__")
        or Path.cwd()
    )
    location = Path(location).resolve()
    start_dir = location if location.is_dir() else location.parent

    for path in [start_dir, *start_dir.parents]:
        if (path / marker).exists():
            return path

    raise RuntimeError(f"Could not find project root above {start_dir}")
