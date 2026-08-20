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
    assert "nucleosuite_version\t0.8.17" in text
    assert "parameter.alpha\t4" in text
    assert "parameter.beta\tvalue" in text
    assert "invocation\tnucleosuite example --alpha 4" in text
