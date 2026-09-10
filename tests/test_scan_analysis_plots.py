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


class _N3Fcn:
    def CptEne(self, ang):
        a = np.asarray(ang, dtype=float)
        return np.cos(np.deg2rad(3.0 * a))


def test_instance_sum_triples_n3_c3_rotor():
    angles = np.linspace(0.0, 350.0, 36)
    v_one = sa.instance_sum_mm_dihed(_N3Fcn(), angles, [0.0])
    v_sum = sa.instance_sum_mm_dihed(_N3Fcn(), angles, [0.0, 120.0, -120.0])
    np.testing.assert_allclose(v_sum, 3.0 * v_one, atol=1e-8)
    assert float(np.max(v_sum) - np.min(v_sum)) == pytest.approx(
        3.0 * float(np.max(v_one) - np.min(v_one)), abs=1e-8
    )


def test_isolate_adds_back_instance_sum_not_one_quartet():
    angles = np.linspace(0.0, 350.0, 36)
    offsets = [0.0, 120.0, -120.0]
    v_sum = sa.instance_sum_mm_dihed(_N3Fcn(), angles, offsets)
    v_one = sa.instance_sum_mm_dihed(_N3Fcn(), angles, [0.0])
    other = 5.0 + 8.0 * np.sin(np.deg2rad(angles))
    e_ll = other + v_sum
    leftover = 0.2 * np.cos(np.deg2rad(angles))
    e_hl = other + leftover
    _a, v_target, v_mm = sa.isolate_dihedral_profiles(
        angles, e_hl, angles, e_ll, v_sum
    )
    np.testing.assert_allclose(v_mm, v_sum, atol=1e-8)
    np.testing.assert_allclose(v_target, leftover, atol=1e-8)
    _a2, v_wrong, _v1 = sa.isolate_dihedral_profiles(
        angles, e_hl, angles, e_ll, v_one
    )
    assert float(np.max(np.abs(v_wrong - leftover))) > 0.5


def test_plot_comparison_extra_curve_and_caption(tmp_path: Path):
    pytest.importorskip("matplotlib")
    angles = np.linspace(0.0, 330.0, 12)
    leftover = 4.0 * (1.0 - np.cos(np.deg2rad(angles)))
    v_sum = 0.6 * (1.0 + np.cos(np.deg2rad(3.0 * angles)))
    v_one = v_sum / 3.0
    cmp = sa.compare_scans(angles, leftover, angles, v_sum)
    out = tmp_path / "compare_xtb_vs_it01_1-2-3-4_dihed.png"
    sa.plot_comparison(
        angles,
        leftover,
        angles,
        v_sum,
        cmp,
        out_path=out,
        hl_label="leftover (HL − MM without type)",
        ll_label="MM DIHE (instance-sum)",
        ylabel="Dihedral term (kcal/mol, min-shifted)",
        extra_curves=[
            (angles, v_one, "one quartet (×3 on this bond)", {"linestyle": "--"})
        ],
        caption="ninst=3  instance-sum ptp=1.20  leftover ptp=8.00",
    )
    assert out.is_file()
    assert out.stat().st_size > 0
