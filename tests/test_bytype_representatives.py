"""Scan only a few bonds per atom-type family when fitting bytype."""

from ffpopt.Workflows import (
    _effective_scans_per_type,
    _select_bytype_representatives,
)


def _rec(a, b, types):
    return {"bond": [a, b], "scan": None, "types": list(types)}


def test_default_two_scans_per_type_when_bytype():
    assert _effective_scans_per_type(True, None) == 2
    assert _effective_scans_per_type(False, None) == 0
    assert _effective_scans_per_type(True, 0) == 0
    assert _effective_scans_per_type(True, 1) == 1


def test_alkyl_tail_collapses_to_two_bonds():
    alkyl = "LIG_c3-c3-c3-c3"
    records = [_rec(i, i + 1, [alkyl, "LIG_c3-c3-c3-hc"]) for i in range(12)]
    selected = _select_bytype_representatives(records, 2)
    assert [r["bond"] for r in selected] == [[0, 1], [1, 2]]


def test_unique_headgroup_type_still_gets_a_scan():
    alkyl = "LIG_c3-c3-c3-c3"
    records = [_rec(i, i + 1, [alkyl]) for i in range(10)]
    records.append(_rec(20, 21, ["LIG_c3-c3-s6-o", alkyl]))
    selected = _select_bytype_representatives(records, 2)
    bonds = [tuple(r["bond"]) for r in selected]
    assert (0, 1) in bonds
    assert (1, 2) in bonds
    assert (20, 21) in bonds
    assert len(selected) == 3


def test_nonpositive_keeps_every_bond():
    records = [_rec(i, i + 1, ["LIG_c3-c3-c3-c3"]) for i in range(5)]
    assert _select_bytype_representatives(records, 0) == records
