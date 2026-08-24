"""Autocomplete, availability counts, and section-title projection."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

from .availability_titles import AvailabilityCounts, title_html


@dataclass
class AvailabilityBindings:
    availability_service: Any
    filtering_pairs: Callable[[], Any]
    original_contigs: list[str]
    original_samples: list[str]
    sample_scope: Any
    widgets: dict[str, Any]
    sort_sample_select: Any
    update_completions: Callable[[Any, list[str]], None]
    total_coverage_count: int
    filtering_title: Any
    contig_title: Any
    sample_title: Any
    mag_title: Any


class AvailabilityController:
    def __init__(self, bindings: AvailabilityBindings) -> None:
        availability_service = bindings.availability_service
        get_filtering_filtered_pairs = bindings.filtering_pairs
        orig_contigs = bindings.original_contigs
        orig_samples = bindings.original_samples
        views = bindings.sample_scope
        widgets = bindings.widgets
        mag_params_sort_sample_select = bindings.sort_sample_select
        update_widget_completions = bindings.update_completions
        _total_coverage_count = bindings.total_coverage_count
        filtering_title = bindings.filtering_title
        contig_title = bindings.contig_title
        sample_title = bindings.sample_title
        mag_title = bindings.mag_title

        def _compute_contig_completions(search_term=""):
            """Core logic — does NOT check the lock. search_term additionally
            narrows results server-side (used for SearchableSelect type-ahead)."""
            search_term = (search_term or "").strip()
            # Tie-break so the currently-selected contig always survives LIMIT
            # 100 (it's a no-op when nothing is selected, since '' never
            # matches a real contig name) — keeps the cap as a real payload
            # guard while never silently dropping what's already selected.
            preserve_contig = widgets["contig_select"].value
            is_mag_view = widgets["has_mags"] and widgets["view_radio"].active == 0
            if is_mag_view:
                # MAG mode: contig list is a child of the selected MAG, not of the sample.
                sel_mag = widgets["mag_select"].value
                if sel_mag and sel_mag in widgets["mag_to_contigs"]:
                    completions = list(widgets["mag_to_contigs"][sel_mag])
                else:
                    completions = list(orig_contigs)
                if search_term:
                    _st = search_term.lower()
                    completions = [c for c in completions if _st in c.lower()]

                # Apply Filtering query builder filters (MAG mode)
                filtered = get_filtering_filtered_pairs()
                if filtered is not None:
                    fsql, fparams = filtered["sql"], filtered["params"]
                    sample_id = None
                    if views.active == 0 and widgets["sample_select"].value:
                        sample_id = widgets["sample_name_to_id"].get(widgets["sample_select"].value)
                    filter_allowed = set(
                        availability_service.filtered_contigs(
                            fsql,
                            fparams,
                            sample_id=sample_id,
                            mag_name=sel_mag if sel_mag in widgets["mag_to_contigs"] else None,
                            search_term=search_term,
                            preserve=preserve_contig,
                        )
                    )
                    completions = [c for c in completions if c in filter_allowed]
            else:
                # Contig mode: filter contigs by selected sample (ONE SAMPLE view).
                if views.active == 0 and widgets["sample_select"].value:
                    sel_sid = widgets["sample_name_to_id"].get(widgets["sample_select"].value)
                    if sel_sid is not None:
                        completions = list(
                            availability_service.contigs_for_sample(sel_sid, search_term, preserve_contig)
                        )
                    else:
                        completions = []
                else:
                    completions = list(orig_contigs)
                    if search_term:
                        _st = search_term.lower()
                        completions = [c for c in completions if _st in c.lower()]

                # Apply Filtering query builder filters (contig mode)
                filtered = get_filtering_filtered_pairs()
                if filtered is not None:
                    fsql, fparams = filtered["sql"], filtered["params"]
                    sample_id = None
                    if views.active == 0 and widgets["sample_select"].value:
                        sample_id = widgets["sample_name_to_id"].get(widgets["sample_select"].value)
                    filter_allowed = set(
                        availability_service.filtered_contigs(
                            fsql,
                            fparams,
                            sample_id=sample_id,
                            search_term=search_term,
                            preserve=preserve_contig,
                        )
                    )
                    completions = [c for c in completions if c in filter_allowed]

            return completions

        def refresh_contig_options_unlocked(search_term=""):
            update_widget_completions(widgets["contig_select"], _compute_contig_completions(search_term))

        def _compute_sample_completions(search_term=""):
            """Core logic — does NOT check the lock. search_term additionally
            narrows results server-side (used for SearchableSelect type-ahead)."""
            search_term = (search_term or "").strip()
            # Tie-break so whatever's currently selected in either consumer of
            # this list (sample_select, and the sort-sample widget which mirrors
            # it) always survives LIMIT 100 — no-op when either is "".
            preserve_sample = widgets["sample_select"].value
            preserve_sort_sample = mag_params_sort_sample_select.value
            is_mag_view = widgets["has_mags"] and widgets["view_radio"].active == 0
            if is_mag_view:
                # MAG mode: filter samples by selected MAG.
                sel_mag = widgets["mag_select"].value
                if views.active == 0 and sel_mag and sel_mag in widgets["mag_to_sample_ids"]:
                    allowed_sids = widgets["mag_to_sample_ids"][sel_mag]
                    _s2id = widgets["sample_name_to_id"]
                    completions = [s for s in orig_samples if _s2id.get(s) in allowed_sids]
                else:
                    completions = list(orig_samples)
                if search_term:
                    _st = search_term.lower()
                    completions = [s for s in completions if _st in s.lower()]

                # Apply Filtering query builder filters (MAG mode)
                filtered = get_filtering_filtered_pairs()
                if filtered is not None:
                    fsql, fparams = filtered["sql"], filtered["params"]
                    contig_id = None
                    if views.active == 0 and widgets["contig_select"].value:
                        contig_id = widgets["contig_name_to_id"].get(widgets["contig_select"].value)
                    filter_allowed = set(
                        availability_service.filtered_samples(
                            fsql,
                            fparams,
                            contig_id=contig_id,
                            mag_name=(
                                sel_mag if contig_id is None and sel_mag in widgets["mag_to_sample_ids"] else None
                            ),
                            search_term=search_term,
                            preserve=(preserve_sample, preserve_sort_sample),
                        )
                    )
                    completions = [s for s in completions if s in filter_allowed]
            else:
                # Contig mode: filter samples by selected contig (ONE SAMPLE view).
                if views.active == 0 and widgets["contig_select"].value:
                    sel_cid = widgets["contig_name_to_id"].get(widgets["contig_select"].value)
                    if sel_cid is not None:
                        completions = list(availability_service.samples_for_contig(sel_cid, search_term))
                    else:
                        completions = []
                else:
                    completions = list(orig_samples)
                    if search_term:
                        _st = search_term.lower()
                        completions = [s for s in completions if _st in s.lower()]

                # Apply Filtering query builder filters (contig mode)
                filtered = get_filtering_filtered_pairs()
                if filtered is not None:
                    fsql, fparams = filtered["sql"], filtered["params"]
                    contig_id = None
                    if views.active == 0 and widgets["contig_select"].value:
                        contig_id = widgets["contig_name_to_id"].get(widgets["contig_select"].value)
                    filter_allowed = set(
                        availability_service.filtered_samples(
                            fsql,
                            fparams,
                            contig_id=contig_id,
                            search_term=search_term,
                            preserve=(preserve_sample, preserve_sort_sample),
                        )
                    )
                    completions = [s for s in completions if s in filter_allowed]
            return completions

        def _sort_sample_completions_for(base_completions, sel_mag=None):
            """Narrow an already-computed sample list to the ones containing the
            currently selected MAG — cheap in-memory filter, no extra query."""
            if not widgets["has_mags"]:
                return base_completions
            sel_mag = sel_mag if sel_mag is not None else widgets["mag_select"].value
            if sel_mag and sel_mag in widgets["mag_to_sample_ids"]:
                allowed_sids = widgets["mag_to_sample_ids"][sel_mag]
                _s2id = widgets["sample_name_to_id"]
                return [s for s in base_completions if _s2id.get(s) in allowed_sids]
            return base_completions

        def refresh_sample_options_unlocked(search_term=""):
            completions = _compute_sample_completions(search_term)
            update_widget_completions(widgets["sample_select"], completions)
            if widgets["has_mags"]:
                update_widget_completions(mag_params_sort_sample_select, _sort_sample_completions_for(completions))

        def _compute_mag_completions(search_term=""):
            """Recompute MAG dropdown based on selected sample (ONE SAMPLE view).
            A MAG is valid iff at least one of its member contigs has coverage for that sample.
            Only meaningful in MAG-mode DBs; no-op otherwise. search_term additionally
            narrows results server-side (used for SearchableSelect type-ahead)."""
            if not widgets["has_mags"]:
                return []
            search_term = (search_term or "").strip()
            mag_to_contigs = widgets["mag_to_contigs"]
            sel_sample = widgets["sample_select"].value
            if views.active == 0 and sel_sample:
                sel_sid = widgets["sample_name_to_id"].get(sel_sample)
                if sel_sid is not None:
                    completions = list(availability_service.mags_for_sample(sel_sid, search_term))
                else:
                    completions = []
            else:
                completions = sorted(mag_to_contigs.keys())
                if search_term:
                    _st = search_term.lower()
                    completions = [m for m in completions if _st in m.lower()]

            # Apply Filtering query builder filters: keep a MAG only if at least
            # one of its member contigs survives the filter.
            filtered = get_filtering_filtered_pairs()
            if filtered is not None:
                fsql, fparams = filtered["sql"], filtered["params"]
                sample_id = None
                if views.active == 0 and sel_sample:
                    sample_id = widgets["sample_name_to_id"].get(sel_sample)
                filter_allowed = set(
                    availability_service.filtered_mags(fsql, fparams, sample_id=sample_id, search_term=search_term)
                )
                completions = [m for m in completions if m in filter_allowed]
            return completions

        def refresh_mag_options_unlocked(search_term=""):
            if not widgets["has_mags"]:
                return
            update_widget_completions(widgets["mag_select"], _compute_mag_completions(search_term))

        _title_fingerprint = {"last": None}

        def invalidate_titles():
            _title_fingerprint["last"] = None

        def _count_contigs_available():
            """True (uncapped) count of contigs valid for the current selection —
            mirrors _compute_contig_completions' base filter (Sample in non-MAG
            mode, MAG-only in MAG mode — MAG mode's Sample->MAG->Contig chain is
            intentional indirection, not filtered here directly), minus the
            LIMIT 100 dropdown-payload cap."""
            if widgets["contig_select"].value:
                return 1
            is_mag_view = widgets["has_mags"] and widgets["view_radio"].active == 0
            if is_mag_view:
                sel_mag = widgets["mag_select"].value
                if sel_mag and sel_mag in widgets["mag_to_contigs"]:
                    return len(widgets["mag_to_contigs"][sel_mag])
                return len(orig_contigs)
            if views.active == 0 and widgets["sample_select"].value:
                sel_sid = widgets["sample_name_to_id"].get(widgets["sample_select"].value)
                if sel_sid is None:
                    return 0
                return availability_service.count_contigs_for_sample(sel_sid)
            return len(orig_contigs)

        def _count_samples_available():
            """Mirror of _count_contigs_available for the Samples header."""
            is_mag_view = widgets["has_mags"] and widgets["view_radio"].active == 0
            if is_mag_view:
                sel_mag = widgets["mag_select"].value
                if views.active == 0 and sel_mag and sel_mag in widgets["mag_to_sample_ids"]:
                    return len(widgets["mag_to_sample_ids"][sel_mag])
                return len(orig_samples)
            if views.active == 0 and widgets["contig_select"].value:
                sel_cid = widgets["contig_name_to_id"].get(widgets["contig_select"].value)
                if sel_cid is None:
                    return 0
                return availability_service.count_samples_for_contig(sel_cid)
            return len(orig_samples)

        def update_section_titles():
            """Update Filtering, Contigs, Samples, and MAGs section titles with current counts."""
            opts_c = widgets["contig_select"].options
            opts_s = widgets["sample_select"].options
            val_c = widgets["contig_select"].value
            val_s = widgets["sample_select"].value
            opts_m = widgets["mag_select"].options if widgets["has_mags"] else ()
            val_m = widgets["mag_select"].value if widgets["has_mags"] else None
            fp = (len(opts_c), val_c, len(opts_s), val_s, len(opts_m), val_m)
            if fp == _title_fingerprint["last"]:
                return
            _title_fingerprint["last"] = fp

            # Use evaluated filter counts when a filter is active; otherwise use preloaded maps.
            filtered = get_filtering_filtered_pairs()
            if filtered is not None:
                presences_count = filtered.get("count_pairs", 0)
                mag_pairs_count = filtered.get("count_mag_pairs", 0)
                fsql_t, fparams_t = filtered["sql"], filtered["params"]
                filtered_mag_counts = (
                    availability_service.filtered_mag_counts(fsql_t, fparams_t, val_m) if val_m else None
                )
                if val_c:
                    contigs_count = 1
                elif val_m:
                    contigs_count = filtered_mag_counts[0]
                else:
                    contigs_count = filtered.get("count_contigs", 0)
                if val_m:
                    samples_count = filtered_mag_counts[1]
                else:
                    samples_count = filtered.get("count_samples", 0)
            else:
                # No active filter — true uncapped counts (the autocomplete
                # lists themselves are capped at LIMIT 100 for dropdown payload
                # size, which would otherwise make this header get stuck
                # reporting "100" regardless of selection)
                contigs_count = _count_contigs_available()
                samples_count = _count_samples_available()
                presences_count = _total_coverage_count
                if widgets["has_mags"]:
                    # Use preloaded mag_to_sample_ids (small: ~2837 MAGs × samples)
                    sel_sid = widgets["sample_name_to_id"].get(val_s) if val_s else None
                    sel_mid = val_m if val_m else None
                    if sel_mid and sel_sid is not None:
                        mag_pairs_count = 1 if sel_sid in widgets["mag_to_sample_ids"].get(sel_mid, set()) else 0
                    elif sel_mid:
                        mag_pairs_count = len(widgets["mag_to_sample_ids"].get(sel_mid, set()))
                    elif sel_sid is not None:
                        mag_pairs_count = sum(1 for sids in widgets["mag_to_sample_ids"].values() if sel_sid in sids)
                    else:
                        mag_pairs_count = sum(len(sids) for sids in widgets["mag_to_sample_ids"].values())

            titles = title_html(
                AvailabilityCounts(
                    contigs_count,
                    samples_count,
                    presences_count,
                    len(set(opts_m) - {""}),
                    mag_pairs_count if widgets["has_mags"] else 0,
                ),
                has_mags=widgets["has_mags"],
            )
            for model, key in (
                (filtering_title, "filtering"),
                (contig_title, "contig"),
                (sample_title, "sample"),
                (mag_title, "mag"),
            ):
                if key in titles and model.text != titles[key]:
                    model.text = titles[key]

        self.compute_contigs = _compute_contig_completions
        self.refresh_contigs = refresh_contig_options_unlocked
        self.compute_samples = _compute_sample_completions
        self.sort_samples_for_mag = _sort_sample_completions_for
        self.refresh_samples = refresh_sample_options_unlocked
        self.compute_mags = _compute_mag_completions
        self.refresh_mags = refresh_mag_options_unlocked
        self.update_titles = update_section_titles
        self.invalidate_titles = invalidate_titles
