"""Regularized Fourier torsion fit (no AmberTools / ASE)."""

from __future__ import annotations

import importlib.util
import os
import sys
import types
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

_FFPOPT = Path(__file__).resolve().parents[1] / "src" / "python" / "lib" / "ffpopt"


def _load_ffpopt_submodules():
    if "ffpopt" not in sys.modules:
        pkg = types.ModuleType("ffpopt")
        pkg.__path__ = [str(_FFPOPT)]
        sys.modules["ffpopt"] = pkg
    loaded = {}
    for name in ("DihedFitRegularize", "Dihedrals"):
        mod_name = f"ffpopt.{name}"
        if mod_name in sys.modules:
            loaded[name] = sys.modules[mod_name]
            continue
        path = _FFPOPT / f"{name}.py"
        spec = importlib.util.spec_from_file_location(mod_name, path)
        mod = importlib.util.module_from_spec(spec)
        sys.modules[mod_name] = mod
        spec.loader.exec_module(mod)
        loaded[name] = mod
    return loaded["DihedFitRegularize"], loaded["Dihedrals"]


reg, dihed = _load_ffpopt_submodules()


def test_tikhonov_beats_cancelling_harmonics():
    tikhonov_svd_solve = reg.tikhonov_svd_solve

    phi = np.deg2rad(np.array([0.0, 2.0, 4.0, 6.0, 8.0]))
    A = np.column_stack([np.cos(n * phi) for n in (1, 2, 3)])
    A = A - np.mean(A, axis=0, keepdims=True)
    y = 0.3 * A[:, 0] + 0.02 * np.array([0.0, 1.0, -1.0, 0.4, -0.3])
    y = y - np.mean(y)
    x_raw = np.linalg.pinv(A, rcond=0.0) @ y
    x_r, info = tikhonov_svd_solve(A, y, lam=0.0, rel_cutoff=1.0e-4)
    assert float(np.linalg.norm(x_raw)) > 10.0 * float(np.linalg.norm(x_r)) + 1.0
    assert info["n_kept"] >= 1


def test_clip_dihed_fcs_caps():
    clip_dihed_fcs = reg.clip_dihed_fcs
    dihed_fc_abs_max = reg.dihed_fc_abs_max

    cap = dihed_fc_abs_max()
    assert cap == 25.0
    out = clip_dihed_fcs([6439.0, -1278.0, 1.4], where="test")
    np.testing.assert_allclose(out, [cap, -cap, 1.4])


def test_fourier_nprim_aic_picks_single_harmonic():
    fit_fourier_nprim = reg.fit_fourier_nprim

    angs = np.linspace(0.0, 350.0, 36)
    y = 1.5 * np.cos(np.deg2rad(angs))
    y = y - np.mean(y)
    with patch.dict(os.environ, {"FFPOPT_DIHED_NPRIM_SELECT": "1"}, clear=False):
        dfcn, x, info = fit_fourier_nprim(angs, y, 3, [0, 1, 2, 3], pname="test-n1")
    assert info["nprim"] == 1
    assert len(dfcn.prims) == 1
    assert float(x[0]) > 0.5


def test_barrier_domain_scales_vphi():
    dense_torsion_ptp = reg.dense_torsion_ptp
    solve_regularized_fcs = reg.solve_regularized_fcs
    GetDihedClasses = dihed.GetDihedClasses

    dfcn = GetDihedClasses(idxs=[0, 1, 2, 3])[1][0]
    A = np.array([[1.0], [-1.0], [0.5], [-0.5]])
    y = np.array([80.0, -80.0, 40.0, -40.0])
    with patch.dict(
        os.environ,
        {
            "FFPOPT_DIHED_BARRIER_ABS": "30",
            "FFPOPT_DIHED_BARRIER_ALPHA": "1",
            "FFPOPT_DIHED_RIDGE_LAMBDA": "0",
        },
        clear=False,
    ):
        x, info = solve_regularized_fcs(A, y, dfcn=dfcn, where="test-barrier")
    assert dense_torsion_ptp(dfcn) <= 30.0 * 1.05
    assert float(np.max(np.abs(x))) <= 25.0
    assert info.get("dense_ptp") is not None


def test_sp3_rotor_policy_zero_and_alkane_cap():
    apply_chemical_rotor_policy = reg.apply_chemical_rotor_policy
    apply_sp3_rotor_policy = reg.apply_sp3_rotor_policy
    classify_dihed_rotor = reg.classify_dihed_rotor
    dense_torsion_ptp = reg.dense_torsion_ptp
    parse_dihed_type_key = reg.parse_dihed_type_key
    GetDihedClasses = dihed.GetDihedClasses

    assert apply_chemical_rotor_policy is apply_sp3_rotor_policy
    assert parse_dihed_type_key("c -ns-c3-c3") == ("c", "ns", "c3", "c3")
    assert classify_dihed_rotor("c3-c3-c3-c3") == "alkane"
    assert classify_dihed_rotor("hc-c3-c3-hc") == "alkane"
    assert classify_dihed_rotor("c3-c3-s6-o") == "sulfate_phosphate"
    assert classify_dihed_rotor("h1-c3-s6-o") == "sulfate_phosphate"
    assert classify_dihed_rotor("c3-c3-n4-c3") == "amine_ammonium"
    assert classify_dihed_rotor("h1-c3-oh-ho") == "alcohol_ether"
    assert classify_dihed_rotor("c3-c3-ss-c3") == "polar_sp3"
    assert classify_dihed_rotor("c3-c3-c3-oh") == "sp3_sp3"
    assert classify_dihed_rotor("o -c -ns-c3") == "unsaturated"
    assert classify_dihed_rotor("h1-c3-c6-c6") == "alkane"
    assert classify_dihed_rotor("h1-c6-c6-h1") == "alkane"
    assert classify_dihed_rotor("oh-c3-c6-os") == "sp3_sp3"
    assert classify_dihed_rotor("oh-c3-c6-h1") == "sp3_sp3"
    assert classify_dihed_rotor("c6-c6-os-c3") == "alcohol_ether"
    assert classify_dihed_rotor("os-c6-os-c3") == "alcohol_ether"
    assert classify_dihed_rotor("h2-c6-os-c3") == "alcohol_ether"
    assert classify_dihed_rotor("h1-c6-os-c6") == "alcohol_ether"
    assert parse_dihed_type_key("CHA_c3-c3-c3-h1") == ("c3", "c3", "c3", "h1")
    assert classify_dihed_rotor("CHA_c3-c3-c3-c3") == "alkane"
    assert classify_dihed_rotor("CHA_c3-c3-s6-o") == "sulfate_phosphate"
    assert classify_dihed_rotor("CHA_c3-c3-n4-c3") == "amine_ammonium"
    assert classify_dihed_rotor("CHA_o-c-ns-c3") == "unsaturated"

    policy_env = {
        "FFPOPT_DIHED_SP3_BARRIER_MAX": "20",
        "FFPOPT_DIHED_ALKANE_BARRIER_MAX": "5",
        "FFPOPT_DIHED_POLAR_SP3_BARRIER_MAX": "8",
        "FFPOPT_DIHED_SULFATE_BARRIER_MAX": "10",
        "FFPOPT_DIHED_SULFATE_BARRIER_CAP": "4",
    }

    stiff = GetDihedClasses(idxs=[0, 1, 2, 3])[1][0]
    stiff.SetFCs([15.4])
    with patch.dict(os.environ, policy_env, clear=False):
        out, action, ptp = apply_sp3_rotor_policy(stiff, "CHA_c3-c3-s6-o", where="test")
    assert action == "cap_sulfate_phosphate"
    assert dense_torsion_ptp(out) <= 4.0 * 1.05
    assert ptp > 10.0

    alk = GetDihedClasses(idxs=[0, 1, 2, 3])[1][0]
    alk.SetFCs([6.0])
    with patch.dict(os.environ, policy_env, clear=False):
        out, action, _ptp = apply_sp3_rotor_policy(alk, "c3-c3-c3-c3", where="test")
    assert action == "cap_alkane"
    assert dense_torsion_ptp(out) <= 5.0 * 1.05

    amine = GetDihedClasses(idxs=[0, 1, 2, 3])[1][0]
    amine.SetFCs([6.0])
    with patch.dict(os.environ, policy_env, clear=False):
        out, action, _ptp = apply_sp3_rotor_policy(amine, "c3-c3-n4-c3", where="test")
    assert action == "cap_amine_ammonium"
    assert dense_torsion_ptp(out) <= 8.0 * 1.05

    sul_cap = GetDihedClasses(idxs=[0, 1, 2, 3])[1][0]
    sul_cap.SetFCs([3.5])
    with patch.dict(os.environ, policy_env, clear=False):
        out, action, _ptp = apply_sp3_rotor_policy(sul_cap, "c3-c3-s6-o", where="test")
    assert action == "cap_sulfate_phosphate"
    assert dense_torsion_ptp(out) <= 4.0 * 1.05

    generic = GetDihedClasses(idxs=[0, 1, 2, 3])[1][0]
    generic.SetFCs([6.0])
    with patch.dict(os.environ, policy_env, clear=False):
        out, action, _ptp = apply_sp3_rotor_policy(
            generic, "c3-c3-c3-oh", where="test"
        )
    assert action == "keep"
    assert float(out.prims[0].fc) == pytest.approx(6.0)

    sugar_h = GetDihedClasses(idxs=[0, 1, 2, 3])[1][0]
    sugar_h.SetFCs([16.3])
    with patch.dict(os.environ, policy_env, clear=False):
        out, action, ptp = apply_sp3_rotor_policy(
            sugar_h, "h1-c3-c6-c6", where="test"
        )
    assert action == "cap_alkane"
    assert dense_torsion_ptp(out) <= 5.0 * 1.05
    assert ptp > 20.0

    sugar_alk = GetDihedClasses(idxs=[0, 1, 2, 3])[1][0]
    sugar_alk.SetFCs([6.0])
    with patch.dict(os.environ, policy_env, clear=False):
        out, action, _ptp = apply_sp3_rotor_policy(
            sugar_alk, "h1-c3-c6-c6", where="test"
        )
    assert action == "cap_alkane"
    assert dense_torsion_ptp(out) <= 5.0 * 1.05

    sugar_cc = GetDihedClasses(idxs=[0, 1, 2, 3])[1][0]
    sugar_cc.SetFCs([15.3])
    with patch.dict(os.environ, policy_env, clear=False):
        out, action, ptp = apply_sp3_rotor_policy(
            sugar_cc, "oh-c3-c6-os", where="test"
        )
    assert action == "cap_sp3_sp3"
    assert dense_torsion_ptp(out) <= 20.0 * 1.05
    assert ptp > 20.0

    sugar_os = GetDihedClasses(idxs=[0, 1, 2, 3])[1][0]
    sugar_os.SetFCs([6.0])
    with patch.dict(os.environ, policy_env, clear=False):
        out, action, _ptp = apply_sp3_rotor_policy(
            sugar_os, "c6-c6-os-c3", where="test"
        )
    assert action == "cap_alcohol_ether"
    assert dense_torsion_ptp(out) <= 8.0 * 1.05

    amide = GetDihedClasses(idxs=[0, 1, 2, 3])[1][0]
    amide.SetFCs([15.0])
    with patch.dict(os.environ, policy_env, clear=False):
        out, action, _ptp = apply_sp3_rotor_policy(amide, "o-c-ns-c3", where="test")
    assert action == "keep"
    assert float(out.prims[0].fc) == pytest.approx(15.0)


def test_scale_fcs_to_leftover_ptp():
    dense_torsion_ptp = reg.dense_torsion_ptp
    scale_fcs_to_ptp = reg.scale_fcs_to_ptp
    GetDihedClasses = dihed.GetDihedClasses

    dfcn = GetDihedClasses(idxs=[0, 1, 2, 3])[1][0]
    dfcn.SetFCs([20.0])
    assert dense_torsion_ptp(dfcn) > 30.0
    scale_fcs_to_ptp(dfcn, 5.0)
    assert dense_torsion_ptp(dfcn) == pytest.approx(5.0, rel=0.02)
