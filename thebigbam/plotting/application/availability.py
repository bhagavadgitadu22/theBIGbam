"""Autocomplete, availability counts, and section-title projection."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

from ..shared.defaults import AUTOCOMPLETE_LIMIT
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


@dataclass(frozen=True)
class AvailabilityResult:
    """One authoritative projection for an autocomplete and its title."""

    options: tuple[str, ...]
    total_count: int


@dataclass(frozen=True)
class AvailabilityScope:
    one_sample: bool
    sample_name: str | None
    sample_id: int | None
    mag_name: str | None
    contig_name: str | None
    contig_id: int | None


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

        _resolution_cache: dict[tuple[Any, ...], AvailabilityResult] = {}

        def _bounded(values, search_term="", preserve=""):
            term = (search_term or "").strip().casefold()
            matches = [value for value in values if not term or term in value.casefold()]
            if preserve in matches:
                matches.remove(preserve)
                matches.insert(0, preserve)
            return matches[:AUTOCOMPLETE_LIMIT]

        def _search(entity: str, fallback, search_term="", preserve=""):
            if not fallback:
                return []
            resolver = getattr(availability_service, entity, None)
            return (
                list(resolver(search_term, preserve))
                if resolver is not None
                else _bounded(fallback, search_term, preserve)
            )

        def _scope() -> AvailabilityScope:
            has_mags = widgets["has_mags"]
            mag_name = widgets["mag_select"].value if has_mags else None
            if not mag_name or mag_name not in widgets.get("mag_to_contigs", {}):
                mag_name = None
            contig_name = widgets["contig_select"].value or None
            if contig_name not in widgets["contig_name_to_id"]:
                contig_name = None
            if mag_name and contig_name not in widgets["mag_to_contigs"].get(mag_name, ()):
                contig_name = None
            sample_name = widgets["sample_select"].value or None
            if sample_name not in widgets["sample_name_to_id"]:
                sample_name = None
            return AvailabilityScope(
                one_sample=views.active == 0,
                sample_name=sample_name,
                sample_id=widgets["sample_name_to_id"].get(sample_name),
                mag_name=mag_name,
                contig_name=contig_name,
                contig_id=widgets["contig_name_to_id"].get(contig_name),
            )

        def _cache_key(entity: str, search_term: str, scope: AvailabilityScope) -> tuple[Any, ...]:
            return (entity, search_term, getattr(availability_service, "revision", 0), scope)

        def _resolve_contigs(search_term=""):
            """Core logic — does NOT check the lock. search_term additionally
            narrows results server-side (used for SearchableSelect type-ahead)."""
            search_term = (search_term or "").strip()
            scope = _scope()
            cache_key = _cache_key("contigs", search_term, scope)
            if cache_key in _resolution_cache:
                return _resolution_cache[cache_key]
            # Tie-break so the currently-selected contig always survives LIMIT
            # 100 (it's a no-op when nothing is selected, since '' never
            # matches a real contig name) — keeps the cap as a real payload
            # guard while never silently dropping what's already selected.
            preserve_contig = widgets["contig_select"].value
            sample_id = scope.sample_id if scope.one_sample else None
            if sample_id is not None:
                completions = list(
                    availability_service.contigs_for_sample(
                        sample_id, search_term, preserve_contig, scope.mag_name
                    )
                )
            elif scope.mag_name:
                completions = _bounded(
                    widgets["mag_to_contigs"][scope.mag_name], search_term, preserve_contig
                )
            else:
                completions = _search("contigs", orig_contigs, search_term, preserve_contig)

            filtered = get_filtering_filtered_pairs()
            if filtered is not None:
                fsql, fparams = filtered["sql"], filtered["params"]
                completions = list(
                    availability_service.filtered_contigs(
                        fsql,
                        fparams,
                        sample_id=sample_id,
                        mag_name=scope.mag_name,
                        search_term=search_term,
                        preserve=preserve_contig,
                    )
                )
                if sample_id is None and scope.mag_name is None:
                    total_count = filtered.get("count_contigs", 0)
                else:
                    total_count = availability_service.count_filtered_contigs(
                        fsql, fparams, sample_id=sample_id, mag_name=scope.mag_name
                    )
            elif sample_id is not None:
                total_count = availability_service.count_contigs_for_sample(sample_id, scope.mag_name)
            elif scope.mag_name:
                total_count = len(widgets["mag_to_contigs"][scope.mag_name])
            else:
                total_count = len(orig_contigs)
            result = AvailabilityResult(tuple(completions), total_count)
            _resolution_cache[cache_key] = result
            return result

        def _compute_contig_completions(search_term=""):
            return list(_resolve_contigs(search_term).options)

        def refresh_contig_options_unlocked(search_term=""):
            update_widget_completions(widgets["contig_select"], _compute_contig_completions(search_term))

        def _resolve_samples(search_term=""):
            """Core logic — does NOT check the lock. search_term additionally
            narrows results server-side (used for SearchableSelect type-ahead)."""
            search_term = (search_term or "").strip()
            scope = _scope()
            cache_key = _cache_key("samples", search_term, scope)
            if cache_key in _resolution_cache:
                return _resolution_cache[cache_key]
            # Tie-break so whatever's currently selected in either consumer of
            # this list (sample_select, and the sort-sample widget which mirrors
            # it) always survives LIMIT 100 — no-op when either is "".
            preserve_sample = widgets["sample_select"].value
            preserve_sort_sample = mag_params_sort_sample_select.value
            contig_id = scope.contig_id if scope.one_sample else None
            scoped_mag = scope.mag_name if contig_id is None else None
            if contig_id is not None:
                completions = list(
                    availability_service.samples_for_contig(contig_id, search_term, preserve_sample)
                )
            elif scope.one_sample and scope.mag_name:
                allowed_sids = widgets["mag_to_sample_ids"][scope.mag_name]
                completions = [
                    sample
                    for sample in orig_samples
                    if widgets["sample_name_to_id"].get(sample) in allowed_sids
                ]
                completions = _bounded(completions, search_term, preserve_sample)
            else:
                completions = _search("samples", orig_samples, search_term, preserve_sample)

            filtered = get_filtering_filtered_pairs()
            if filtered is not None:
                fsql, fparams = filtered["sql"], filtered["params"]
                # Applied MAG filters remain scoped to a selected MAG in ALL
                # mode as well, even though the Samples selector is hidden.
                if contig_id is None:
                    scoped_mag = scope.mag_name
                completions = list(
                    availability_service.filtered_samples(
                        fsql,
                        fparams,
                        contig_id=contig_id,
                        mag_name=scoped_mag,
                        search_term=search_term,
                        preserve=(preserve_sample, preserve_sort_sample),
                    )
                )
                if contig_id is None and scoped_mag is None:
                    total_count = filtered.get("count_samples", 0)
                else:
                    total_count = availability_service.count_filtered_samples(
                        fsql, fparams, contig_id=contig_id, mag_name=scoped_mag
                    )
            elif contig_id is not None:
                total_count = availability_service.count_samples_for_contig(contig_id)
            elif scope.one_sample and scope.mag_name:
                total_count = len(widgets["mag_to_sample_ids"][scope.mag_name])
            else:
                total_count = len(orig_samples)
            result = AvailabilityResult(tuple(completions), total_count)
            _resolution_cache[cache_key] = result
            return result

        def _compute_sample_completions(search_term=""):
            return list(_resolve_samples(search_term).options)

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

        def _resolve_mags(search_term=""):
            """Recompute MAG dropdown based on selected sample (ONE SAMPLE view).
            A MAG is valid iff at least one of its member contigs has coverage for that sample.
            Only meaningful in MAG-mode DBs; no-op otherwise. search_term additionally
            narrows results server-side (used for SearchableSelect type-ahead)."""
            if not widgets["has_mags"]:
                return AvailabilityResult((), 0)
            search_term = (search_term or "").strip()
            scope = _scope()
            cache_key = _cache_key("mags", search_term, scope)
            if cache_key in _resolution_cache:
                return _resolution_cache[cache_key]
            mag_to_contigs = widgets["mag_to_contigs"]
            preserve_mag = widgets["mag_select"].value
            if scope.one_sample and scope.sample_id is not None:
                completions = list(
                    availability_service.mags_for_sample(scope.sample_id, search_term, preserve_mag)
                )
            else:
                completions = _search("mags", mag_to_contigs, search_term, preserve_mag)

            # Apply Filtering query builder filters: keep a MAG only if at least
            # one of its member contigs survives the filter.
            filtered = get_filtering_filtered_pairs()
            if filtered is not None:
                fsql, fparams = filtered["sql"], filtered["params"]
                sample_id = scope.sample_id if scope.one_sample else None
                completions = list(
                    availability_service.filtered_mags(
                        fsql,
                        fparams,
                        sample_id=sample_id,
                        search_term=search_term,
                        preserve=preserve_mag,
                    )
                )
                completions = completions[:AUTOCOMPLETE_LIMIT]
            if filtered is not None:
                total_count = availability_service.count_filtered_mags(fsql, fparams, sample_id=sample_id)
            elif scope.one_sample and scope.sample_id is not None:
                # MAG options are capped, so count the full set from the small
                # preloaded MAG-to-sample relation.
                total_count = sum(
                    1
                    for sample_ids in widgets["mag_to_sample_ids"].values()
                    if scope.sample_id in sample_ids
                )
            else:
                total_count = len(mag_to_contigs)
            result = AvailabilityResult(tuple(completions), total_count)
            _resolution_cache[cache_key] = result
            return result

        def _compute_mag_completions(search_term=""):
            return list(_resolve_mags(search_term).options)

        def refresh_mag_options_unlocked(search_term=""):
            if not widgets["has_mags"]:
                return
            update_widget_completions(widgets["mag_select"], _compute_mag_completions(search_term))

        _title_fingerprint = {"last": None}

        def invalidate_titles():
            _title_fingerprint["last"] = None
            _resolution_cache.clear()

        def update_section_titles():
            """Update Filtering, Contigs, Samples, and MAGs section titles with current counts."""
            scope = _scope()
            fp = (scope, getattr(availability_service, "revision", 0))
            if fp == _title_fingerprint["last"]:
                return
            _title_fingerprint["last"] = fp

            contig_result = _resolve_contigs()
            sample_result = _resolve_samples()
            mag_result = _resolve_mags()
            contigs_count = contig_result.total_count
            samples_count = sample_result.total_count

            filtered = get_filtering_filtered_pairs()
            if filtered is not None:
                presences_count = filtered.get("count_pairs", 0)
                mag_pairs_count = filtered.get("count_mag_pairs", 0)
            else:
                presences_count = _total_coverage_count
                if widgets["has_mags"]:
                    # Use preloaded mag_to_sample_ids (small: ~2837 MAGs × samples)
                    sel_sid = scope.sample_id
                    sel_mid = scope.mag_name
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
                    mag_result.total_count,
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
