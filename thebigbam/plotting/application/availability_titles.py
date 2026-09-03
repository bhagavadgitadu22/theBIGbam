"""Pure title formatting for availability counts."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class AvailabilityCounts:
    contigs: int
    samples: int
    presences: int
    mags: int = 0
    mag_pairs: int = 0


def title_html(counts: AvailabilityCounts, *, has_mags: bool) -> dict[str, str]:
    pair_text = f"{counts.presences} contig/sample pairs"
    if has_mags:
        pair_text += f", {counts.mag_pairs} MAG/sample pairs"
    result = {
        "filtering": (
            "<span style='font-size: 1.2em;'><b>Filtering</b></span>"
            f" <span style='font-size: 0.85em;'>{pair_text}</span>"
        ),
        "contig": (
            "<span style='font-size: 1.2em;'><b>Contigs</b></span> "
            f"<span style='font-size: 0.85em;'>{counts.contigs} available</span>"
        ),
        "sample": (
            "<span style='font-size: 1.2em;'><b>Samples</b></span> "
            f"<span style='font-size: 0.85em;'>{counts.samples} available</span>"
        ),
    }
    if has_mags:
        result["mag"] = (
            "<span style='font-size: 1.2em;'><b>MAGs</b></span> "
            f"<span style='font-size: 0.85em;'>{counts.mags} available</span>"
        )
    return result
