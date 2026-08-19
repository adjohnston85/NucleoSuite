from pathlib import Path

from nucleosuite.nrl import default_output_prefix
from nucleosuite.output_naming import compact_parameter, parameterized_prefix


def test_nrl_explicit_prefix_is_a_parameterized_base() -> None:
    assert default_output_prefix(
        Path("sample.tsv"),
        2,
        140,
        1,
        base="Gaffney32_145_147_NRL",
    ) == Path("Gaffney32_145_147_NRL_peakres1_min2_max140")


def test_parameter_tokens_are_safe_and_idempotent() -> None:
    prefix = parameterized_prefix(
        "sample",
        (("match", "many-to-one"), ("maxdist", None), ("lower", -80)),
    )
    assert prefix == Path("sample_matchmany-to-one_maxdistnone_lowerneg80")
    assert parameterized_prefix(
        prefix,
        (("match", "many-to-one"), ("maxdist", None), ("lower", -80)),
    ) == prefix
    assert compact_parameter(1.5) == "1p5"
