from pathlib import Path

import numpy as np

from nucleosuite.peak_fdr import annotate_peak_fdr, empirical_peak_pvalues, empirical_peak_qvalues


def _write_bed(path: Path, scores):
    path.write_text(
        "".join(
            f"chr1\t{index * 10}\t{index * 10 + 8}\tpeak{index}\t{score}\t.\t"
            f"{index * 10 + 4}\t{index * 10 + 5}\n"
            for index, score in enumerate(scores)
        ),
        encoding="utf-8",
    )


def test_empirical_peak_qvalues_use_random_peak_counts_and_are_monotonic():
    observed = np.asarray([30.0, 24.0, 18.0, 12.0])
    randomized = [np.asarray([20.0, 10.0])]
    qvalues, observed_counts, mean_random_counts = empirical_peak_qvalues(
        observed, randomized
    )

    assert observed_counts.tolist() == [1, 2, 3, 4]
    assert mean_random_counts.tolist() == [0.0, 0.0, 1.0, 1.0]
    assert np.all((qvalues >= 0) & (qvalues <= 1))
    assert qvalues[0] <= qvalues[1] <= qvalues[2] <= qvalues[3]
    assert qvalues.tolist() == [0.5, 0.5, 0.5, 0.5]

def test_empirical_peak_pvalues_use_pooled_random_upper_tail():
    observed = np.asarray([30.0, 18.0, 5.0])
    randomized = [np.asarray([20.0, 10.0]), np.asarray([15.0])]
    pvalues = empirical_peak_pvalues(observed, randomized)
    assert np.allclose(pvalues, [0.25, 0.5, 1.0])



def test_peak_fdr_preserves_input_fields_and_appends_pvalue_and_fdr(tmp_path: Path):
    sample = tmp_path / "sample.bed"
    random = tmp_path / "random.bed"
    _write_bed(sample, [30, 24, 18, 12])
    _write_bed(random, [20, 10])

    result = annotate_peak_fdr(
        sample,
        [random],
        fdr_threshold=0.6,
        output_prefix=tmp_path / "result",
    )

    source_fields = [line.split("\t") for line in sample.read_text().splitlines()]
    output_fields = [
        line.split("\t") for line in result.annotated_path.read_text().splitlines()
    ]
    assert len(output_fields) == len(source_fields)
    assert all(output[:-2] == source for output, source in zip(output_fields, source_fields))
    assert all(0 <= float(output[-2]) <= 1 for output in output_fields)
    assert all(0 <= float(output[-1]) <= 1 for output in output_fields)
    assert result.significant_path is not None
    assert result.summary_path.is_file()


def test_peak_fdr_without_cutoff_writes_every_peak_and_no_filtered_bed(tmp_path: Path):
    sample = tmp_path / "sample.bed"
    random = tmp_path / "random.bed"
    _write_bed(sample, [10, 5])
    _write_bed(random, [4])

    result = annotate_peak_fdr(sample, [random])

    assert len(result.annotated_path.read_text().splitlines()) == 2
    assert result.significant_path is None
    assert result.significant_peaks is None


def test_peak_fdr_accepts_an_empty_randomized_peak_callset(tmp_path: Path):
    sample = tmp_path / "sample.bed"
    random = tmp_path / "random.bed"
    _write_bed(sample, [10, 5])
    random.write_text("", encoding="utf-8")

    result = annotate_peak_fdr(sample, [random])

    assert len(result.annotated_path.read_text().splitlines()) == 2
    assert all(0 <= float(line.split("\t")[-1]) <= 1 for line in result.annotated_path.read_text().splitlines())


def test_peak_fdr_accepts_an_empty_observed_peak_callset(tmp_path: Path):
    sample = tmp_path / "sample.bed"
    random = tmp_path / "random.bed"
    sample.write_text("", encoding="utf-8")
    _write_bed(random, [4])

    result = annotate_peak_fdr(sample, [random], fdr_threshold=0.05)

    assert result.annotated_path.read_text() == ""
    assert result.significant_path is not None
    assert result.significant_path.read_text() == ""
