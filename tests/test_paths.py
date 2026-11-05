"""The store root resolves to the main checkout, which git itself names."""

import subprocess
from pathlib import Path

import pytest

from common._paths import MAIN_CHECKOUT, REPO_ROOT


def test_main_checkout_matches_git():
    """MAIN_CHECKOUT is the parent of git's common dir, from a worktree or the main checkout."""
    if not (REPO_ROOT / ".git").exists():
        pytest.skip("not a git checkout")
    git_common_dir = subprocess.run(
        ["git", "rev-parse", "--path-format=absolute", "--git-common-dir"],
        cwd=REPO_ROOT, capture_output=True, text=True, check=True,
    ).stdout.strip()
    expected_checkout = Path(git_common_dir).resolve().parent
    assert MAIN_CHECKOUT == expected_checkout
