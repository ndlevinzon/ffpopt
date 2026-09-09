"""Total-energy vs isolated-DIHE comparison plots (no AmberTools / ASE)."""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest

_SCAN = (
    Path(__file__).resolve().parents[1]
    / "src"
    / "python"
    / "lib"
    / "ffpopt"
    / "ScanAnalysis.py"
)


def _load_scananalysis():
    spec = importlib.util.spec_from_file_location("scananalysis_under_test", _SCAN)
    mod = importlib.util.module_from_spec(spec)
    import sys

    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


sa = _load_scananalysis()


def _v_mm(ang: float) -> float:
    return 2.0 * (1.0 + math.cos(math.radians(ang)))


def test_isolate_recovers_dihedral_when_other_mm_terms_match():
    angles = np.linspace(0.0, 330.0, 12)
    v = np.array([_v_mm(a) for a in angles])
    other = 3.0 + 0.1 * np.sin(np.deg2rad(2.0 * angles))
    e_ll = other + v
    e_hl = other + v
    a, v_target, v_mm = sa.isolate_dihedral_profiles(
        angles, e_hl, angles, e_ll, _v_mm
    )
    np.testing.assert_allclose(a, angles)
    np.testing.assert_allclose(v_target, v_mm, atol=1e-8)
    np.testing.assert_allclose(v_mm, v, atol=1e-8)


def test_isolate_qm_leftover_is_hl_minus_mm_other():
    angles = np.linspace(0.0, 330.0, 12)
    v = np.array([_v_mm(a) for a in angles])
    other = np.full_like(angles, 1.5)
    e_ll = other + v
    extra = 0.4 * np.cos(np.deg2rad(3.0 * angles))
    e_hl = other + extra
    _a, v_target, v_mm = sa.isolate_dihedral_profiles(
        angles, e_hl, angles, e_ll, _v_mm
    )
    np.testing.assert_allclose(v_target, extra, atol=1e-8)
    np.testing.assert_allclose(v_mm, v, atol=1e-8)


def test_plot_comparison_ylabel_total_energy(tmp_path: Path):
    pytest.importorskip("matplotlib")
    angles = np.linspace(0.0, 330.0, 12)
    e_hl = 0.5 * (1.0 - np.cos(np.deg2rad(angles)))
    e_ll = 0.4 * (1.0 - np.cos(np.deg2rad(angles)))
    cmp = sa.compare_scans(angles, e_hl, angles, e_ll)
    out = tmp_path / "compare_xtb_vs_it01_0-1-2-3.png"
    sa.plot_comparison(
        angles,
        e_hl,
        angles,
        e_ll,
        cmp,
        out_path=out,
        hl_label="xtb",
        ll_label="it01",
        ylabel="Total energy (kcal/mol, min-shifted)",
    )
    assert out.is_file()
    assert out.stat().st_size > 0
