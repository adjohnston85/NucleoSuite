from __future__ import annotations

import random
from pathlib import Path


def test_blacklist_merges_aliases_and_filters_fragment_source(tmp_path):
    from nucleosuite.core.blacklist import load_blacklist
    from nucleosuite.core.fragment_inputs import BlacklistFragmentSource, FragmentSource

    bed = tmp_path / "blacklist.bed"
    bed.write_text("1\t10\t20\nchr1\t18\t30\n", encoding="utf-8")
    index = load_blacklist(bed, ["chr1"], [100])
    assert index is not None
    assert index.intervals["chr1"] == ((10, 30),)
    assert index.summary.interval_count == 1
    assert index.summary.blacklisted_bases == 20

    class FakeSource(FragmentSource):
        references = ["chr1"]
        lengths = [100]
        input_paths = ["fake.bed"]
        kind = "fragments"

        def fetch(self, contig, start, end, **kwargs):
            return [(0, 10), (9, 11), (30, 40)]

    source = BlacklistFragmentSource(FakeSource(), index)
    assert source.fetch(
        "chr1", 0, 100, max_per_coordinate=0, subsample=None,
        dedup_scope="all_bams",
    ) == [(0, 10), (30, 40)]
    assert source.fragments_excluded == 1


def test_dinucleotide_randomization_uses_available_other_anchor_and_blacklist(tmp_path):
    from nucleosuite.core.blacklist import load_blacklist
    from nucleosuite.core.randomization import RandomizationBlock, place_dinucleotide_matched

    bed = tmp_path / "blacklist.bed"
    bed.write_text("chr1\t4\t8\n", encoding="utf-8")
    blacklist = load_blacklist(bed, ["chr1"], [24])
    sequence = "AACCAACCAACCAACCAACCAACC"
    block = RandomizationBlock("chr1", 0, len(sequence), sequence, blacklist=blacklist)
    output, status, selected, matched, reason = place_dinucleotide_matched(
        (0, 4),
        start_dinuc="NN",  # deliberately unavailable
        end_dinuc="CC",
        block=block,
        rng=random.Random(4),
        anchor_prob_start=1.0,
        fallback="skip",
    )
    assert selected == "start"
    assert matched == "end"
    assert status == "matched"
    assert reason is None
    assert output is not None and output != (0, 4)
    assert not blacklist.overlaps("chr1", *output)


def test_randomization_rejects_candidate_fragments_spanning_non_acgt():
    from nucleosuite.core.randomization import RandomizationBlock

    block = RandomizationBlock("chr1", 0, 12, "AAAANNAAAAAA")
    candidates = block.anchor_candidates("AA", 4, "start", (0, 4))
    assert candidates
    assert all("N" not in block.sequence[start : start + 4] for start in candidates)
    assert block.non_acgt_candidate_rejections > 0


def test_category_colours_do_not_repeat_for_sixteen_dinucleotides():
    from nucleosuite.plotting import category_colors

    colours = category_colors(16)
    assert len(colours) == 16
    assert len({tuple(colour) for colour in colours}) == 16


def test_blacklist_ignores_contigs_absent_from_split_source(tmp_path):
    from nucleosuite.core.blacklist import load_blacklist

    bed = tmp_path / "whole_assembly_blacklist.bed"
    bed.write_text("chr1\t10\t20\nchr17\t30\t40\n", encoding="utf-8")
    index = load_blacklist(bed, ["chr17"], [100])
    assert index is not None
    assert index.intervals == {"chr17": ((30, 40),)}
    assert index.summary.ignored_interval_count == 1
