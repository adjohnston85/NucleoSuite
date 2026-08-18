from pathlib import Path


def test_randomize_fragments_is_block_bounded_and_transactional():
    text = Path("src/nucleosuite/randomize_fragments_command.py").read_text()
    assert "sqlite3" not in text
    assert "while core_start < contig_length" in text
    assert "randomized_chunk.sort()" in text
    assert "_validate_randomized_bed(" in text
    assert "os.replace(temporary_paths[path], path)" in text
    assert text.index("plot_count_profile(") < text.index("os.replace(temporary_paths[path], path)")


def test_candidate_indexes_do_not_copy_every_candidate_per_fragment():
    text = Path("src/nucleosuite/core/randomization.py").read_text()
    assert 'values = array("Q")' in text
    assert "_uniform_cache" in text
    assert "_anchor_cache" in text
    assert "allowed = len(cached) - int(original_present)" in text


def test_original_position_is_forbidden_in_all_placement_paths():
    from nucleosuite.core.randomization import RandomizationBlock, uniform_randomize_fragment
    import random

    block = RandomizationBlock("chr1", 0, 100, "A" * 100)
    original = (10, 20)
    outputs = [uniform_randomize_fragment(original, block, random.Random(seed)) for seed in range(25)]
    assert all(output is not None and output != original for output in outputs)
