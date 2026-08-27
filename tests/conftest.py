"""Shared test fixtures.

`abnumber_tool.py` and `anarci_tool.py` do module-level `from abnumber import ...`
and `from anarci import anarci`, which are conda/bioconda-only packages not
installed in this Tier 1 environment. We stub minimal fake modules into
sys.modules so those imports succeed at collection time, without ever
exercising real ANARCI/AbNumber behavior (the functions that actually call
into them are Tier 2 scope and are not tested here).
"""

import sys
from pathlib import Path
from types import ModuleType

import pandas as pd
import pytest

THERASABDAB_CSV = Path(__file__).parent.parent / "test_files" / "therasabdab_sample.csv"

def _really_importable(module_name: str) -> bool:
    """True if `module_name` genuinely resolves to a real installed package
    (not a stub this file may have previously inserted into sys.modules)."""
    import importlib

    try:
        importlib.import_module(module_name)
        return not getattr(sys.modules[module_name], "__anarci_toolz_stub__", False)
    except ImportError:
        return False


# Recorded once, *before* any stubbing below, so Tier 2 test files can tell
# a real install apart from our stub via `require_real_module` — a plain
# `pytest.importorskip("abnumber")` would be fooled by the stub, since by
# the time a test module runs, sys.modules["abnumber"] already resolves
# (this bit us: Tier 2 tests silently ran against the stub instead of
# skipping, and a real multiprocessing Pool trying to call/pickle the stub
# lambda hung the whole suite).
ABNUMBER_AVAILABLE = _really_importable("abnumber")
ANARCI_AVAILABLE = _really_importable("anarci")


def require_real_module(module_name: str) -> None:
    """Tier 2 (test_*_bio.py) equivalent of `pytest.importorskip`, aware of
    the Tier 1 stubs below."""
    available = {"abnumber": ABNUMBER_AVAILABLE, "anarci": ANARCI_AVAILABLE}[module_name]
    if not available:
        pytest.skip(f"{module_name} is not installed", allow_module_level=True)


# `abnumber_tool.py` and `anarci_tool.py` do module-level `from abnumber import ...`
# and `from anarci import anarci`, which are conda/bioconda-only packages not
# installed in the Tier 1 environment. We stub minimal fake modules into
# sys.modules so those imports succeed at collection time for Tier 1, without
# ever exercising real ANARCI/AbNumber behavior.
if not ABNUMBER_AVAILABLE:
    fake_abnumber = ModuleType("abnumber")
    fake_abnumber.__anarci_toolz_stub__ = True
    fake_abnumber.Chain = object
    fake_abnumber.Position = object
    sys.modules["abnumber"] = fake_abnumber

if not ANARCI_AVAILABLE:
    fake_anarci = ModuleType("anarci")
    fake_anarci.__anarci_toolz_stub__ = True
    fake_anarci.anarci = lambda *args, **kwargs: (None, None, None)
    sys.modules["anarci"] = fake_anarci


@pytest.fixture
def small_seq_df():
    return pd.DataFrame(
        {
            "sequence_aa": ["EVQLVESGGGLVQPGG", "DIQMTQSPSSLSASVG"],
            "name": ["ab1", "ab2"],
        }
    )


@pytest.fixture(scope="session")
def therasabdab_full_df():
    """The README's own 100-row worked example. Loaded once per session since
    Tier 2 tests are slow (real ANARCI/AbNumber calls)."""
    return pd.read_csv(THERASABDAB_CSV)


@pytest.fixture(scope="session")
def therasabdab_subset_df(therasabdab_full_df):
    """First 8 rows of the fixture, which includes the Timigutuzumab and
    Inebilizumab rows used for README-value spot checks, kept small to bound
    real ANARCI/AbNumber runtime in per-function Tier 2 tests."""
    return therasabdab_full_df.head(8).reset_index(drop=True)
