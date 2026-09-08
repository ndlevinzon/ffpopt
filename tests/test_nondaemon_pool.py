"""Scheduling helpers used by the performance-oriented spawn pools."""

from ffpopt.NondaemonPool import (
    acquire_wavefront_pool,
    close_reused_wavefront_pool,
    split_nproc_for_items,
)


def test_split_prefer_depth_allows_two_d():
    assert split_nproc_for_items(
        32, 4, prefer_depth=True, flatten_nested=False
    ) == (4, 8)


def test_split_nested_flattens_to_one_axis():
    assert split_nproc_for_items(
        32, 4, prefer_depth=True, flatten_nested=True
    ) == (1, 32)


def test_split_single_item_keeps_all_cores():
    assert split_nproc_for_items(16, 1, prefer_depth=True, flatten_nested=False) == (
        1,
        16,
    )


def test_wavefront_pool_reuse_same_nproc():
    pool_a, owns_a = acquire_wavefront_pool(2)
    pool_b, owns_b = acquire_wavefront_pool(2)
    try:
        assert pool_a is pool_b
        assert owns_a is False
        assert owns_b is False
    finally:
        close_reused_wavefront_pool()


def test_clash_precheck_skips_bonded_pairs():
    from ffpopt.WavefrontIpc import has_nonbonded_clash

    pos = [[0.0, 0.0, 0.0], [0.4, 0.0, 0.0], [3.0, 0.0, 0.0]]
    clashed, i, j, dist = has_nonbonded_clash(pos, bonds=[(0, 1)], min_dist=0.8)
    assert clashed is False


def test_clash_precheck_flags_nonbonded():
    from ffpopt.WavefrontIpc import has_nonbonded_clash

    pos = [[0.0, 0.0, 0.0], [0.4, 0.0, 0.0], [3.0, 0.0, 0.0]]
    clashed, i, j, dist = has_nonbonded_clash(pos, bonds=[], min_dist=0.8)
    assert clashed is True
    assert {i, j} == {0, 1}
    assert dist < 0.8


def test_checkpoint_every_env(monkeypatch):
    from ffpopt.WavefrontIpc import wf_checkpoint_every

    monkeypatch.delenv("FFPOPT_WF_CHECKPOINT_EVERY", raising=False)
    assert wf_checkpoint_every(8) == 8
    monkeypatch.setenv("FFPOPT_WF_CHECKPOINT_EVERY", "24")
    assert wf_checkpoint_every(8) == 24


def test_require_main_guard_is_quiet_outside_spawn():
    from ffpopt.WavefrontIpc import require_main_guard_for_spawn

    require_main_guard_for_spawn("test")
