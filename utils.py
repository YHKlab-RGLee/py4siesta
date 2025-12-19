import os
import shutil
from contextlib import contextmanager
from pathlib import Path
from typing import Optional


@contextmanager
def working_dir(path: Path):
    """Temporarily change the working directory to ``path``.

    The previous working directory is restored even if an exception occurs.
    """

    previous_cwd = Path.cwd()
    target = Path(path)
    try:
        target.mkdir(parents=True, exist_ok=True)
        os.chdir(target)
        yield target
    finally:
        os.chdir(previous_cwd)


def copy_contents(src: Path, dst: Path):
    """Copy the contents of ``src`` into ``dst``."""

    src = Path(src)
    dst = Path(dst)
    dst.mkdir(parents=True, exist_ok=True)
    for item in src.iterdir():
        target = dst / item.name
        if item.is_dir():
            shutil.copytree(item, target)
        else:
            shutil.copy2(item, target)


def last_matching_line(path: Path, keyword: str) -> Optional[str]:
    """Return the last line in ``path`` containing ``keyword`` or ``None``."""

    path = Path(path)
    if not path.is_file():
        return None

    with path.open("r", errors="ignore") as file:
        matching = [line for line in file if keyword in line]
    return matching[-1] if matching else None
