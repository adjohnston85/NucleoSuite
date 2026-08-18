"""Shared input/output helpers used throughout NucleoSuite."""

from __future__ import annotations

from pathlib import Path

from .intervals import (
    BIGBED_SUFFIXES,
    INTERVAL_FORMATS,
    convert_bed_to_bigbed,
    finalise_interval_files,
    is_bigbed_path,
    open_interval_text,
)

open_text = open_interval_text


def strip_known_suffix(path: str | Path) -> str:
    """Return a basename with common genomics file suffixes removed."""
    name = Path(path).name
    for suffix in (
        ".bed.gz", ".BED.gz", ".bigBed", ".bigbed", ".bb", ".bed", ".BED",
        ".bigWig", ".bigwig", ".bw",
        ".txt.gz", ".tsv.gz", ".txt", ".tsv", ".gz",
    ):
        if name.endswith(suffix):
            return name[:-len(suffix)]
    return Path(name).stem


__all__ = [
    "BIGBED_SUFFIXES",
    "INTERVAL_FORMATS",
    "convert_bed_to_bigbed",
    "finalise_interval_files",
    "is_bigbed_path",
    "open_text",
    "strip_known_suffix",
]
