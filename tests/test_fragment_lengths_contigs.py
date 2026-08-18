from nucleosuite.fragment_lengths import build_parser, parse_contigs


def test_fragment_lengths_accepts_space_separated_contigs():
    args = build_parser().parse_args(
        ["--fragments", "fragments.bed", "--contigs", "19", "20", "21", "22"]
    )
    assert args.contigs == ["19", "20", "21", "22"]
    assert parse_contigs(args.contigs) == {"19", "20", "21", "22"}


def test_fragment_lengths_accepts_comma_separated_contigs():
    args = build_parser().parse_args(
        ["--fragments", "fragments.bed", "--contigs", "19,20,21,22"]
    )
    assert parse_contigs(args.contigs) == {"19", "20", "21", "22"}


def test_fragment_lengths_all_selector_wins():
    assert parse_contigs(["19", "all", "22"]) is None
