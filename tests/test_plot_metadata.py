from pathlib import Path


def test_save_figure_writes_parameter_metadata(tmp_path: Path) -> None:
    import matplotlib.pyplot as plt
    from nucleosuite.plotting import configure_plot_metadata, save_figure

    configure_plot_metadata(
        "example",
        ["example", "--alpha", "4"],
        {"alpha": 4, "beta": "value", "items": [1, 2]},
    )
    fig, ax = plt.subplots()
    ax.plot([0, 1], [0, 1])
    output = save_figure(fig, tmp_path / "plot.png")
    plt.close(fig)
    metadata = output.with_name(output.stem + "_metadata.tsv")
    text = metadata.read_text()
    assert "nucleosuite_version\t0.9.1" in text
    assert "parameter.alpha\t4" in text
    assert "parameter.beta\tvalue" in text
    assert "invocation\tnucleosuite example --alpha 4" in text


def test_shared_plot_source_metadata_retains_primary_recipe_and_all_associations(tmp_path: Path) -> None:
    from nucleosuite.plotting import configure_plot_metadata, plot_source_metadata_path
    from nucleosuite.profile_plots import plot_dinucleotide_profile, plot_ww_ss_profile

    source = tmp_path / "sample_dinuc_profile.tsv"
    source.write_text(
        "position\tn_valid\tAA_pct\tAT_pct\tWW_pct\tSS_pct\n"
        "-1\t10\t5\t6\t11\t8\n"
        "0\t10\t7\t8\t15\t9\n"
        "1\t10\t5\t6\t11\t8\n",
        encoding="utf-8",
    )
    configure_plot_metadata("dinuc-profile", ["dinuc-profile"], {"major_grid_color": "0.65"})
    first = plot_dinucleotide_profile(str(source), str(tmp_path / "sample_dinuc_profile.png"))
    second = plot_ww_ss_profile(str(source), str(tmp_path / "sample_dinuc_profile_ww_ss.png"))
    assert first is not None and second is not None

    metadata = plot_source_metadata_path(source).read_text(encoding="utf-8")
    assert "detected_plot_type\tdinucleotide-profile" in metadata
    assert "associated_plot_types\tdinucleotide-profile,ww-ss-profile" in metadata
    assert "resolved_title\tsample — Dinucleotide profile" in metadata
    assert "plot.dinucleotide-profile.resolved_title\tsample — Dinucleotide profile" in metadata
    assert "plot.ww-ss-profile.resolved_title\tsample — WW/SS dinucleotide profile" in metadata
