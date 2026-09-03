"""Sample, contig, MAG, search, and genomic-position interaction callbacks."""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from typing import Any, Callable

from bokeh.io import curdoc

from ..controls.searchable_select import decode_search_request

CONTIG_TO_MAG_SAMPLE_ORDER_CATEGORY = {
    "Coverage": "MAG coverage",
    "Misassembly": "MAG misassembly",
    "Microdiversity": "MAG microdiversity",
}
def sample_order_category_for_view(category: str, *, mag_view: bool) -> str:
    """Map equivalent sample-order categories to the active subject view."""
    if mag_view:
        return CONTIG_TO_MAG_SAMPLE_ORDER_CATEGORY.get(category, category)
    return category


def translate_mag_window_to_contig(
    start: int, end: int, *, offset: int, contig_length: int
) -> tuple[int, int]:
    """Translate a 1-based MAG window and keep it inside one contig."""
    if contig_length <= 0 or start < 1 or end < start:
        return 1, max(0, contig_length)
    width = end - start + 1
    if width >= contig_length:
        return 1, contig_length
    local_start = start - offset
    local_end = end - offset
    if local_end < 1 or local_start > contig_length:
        return 1, contig_length
    if local_start < 1:
        local_end += 1 - local_start
        local_start = 1
    if local_end > contig_length:
        local_start -= local_end - contig_length
        local_end = contig_length
    return max(1, local_start), min(contig_length, local_end)


@dataclass(frozen=True)
class SubjectBindings:
    widgets: dict[str, Any]
    interaction_lock: dict[str, bool]
    compute_contigs: Callable[[str], list[str]]
    compute_samples: Callable[[str], list[str]]
    compute_mags: Callable[[str], list[str]]
    push_completions: Callable[[Any, list[str]], None]
    refresh_contigs: Callable[..., None]
    refresh_samples: Callable[..., None]
    refresh_mags: Callable[..., None]
    update_titles: Callable[[], None]
    from_position: Any
    to_position: Any
    sample_order_category: Any
    sample_contig_categories: list[str]
    sample_mag_categories: list[str]
    sample_current_categories: list[str]
    schedule_transition: Callable[[Callable[[], None]], None] = lambda callback: callback()


class SubjectController:
    def __init__(self, bindings: SubjectBindings) -> None:
        self.bindings = bindings

    def attach(self) -> None:
        widgets = self.bindings.widgets
        widgets["mag_select"].param.watch(self.mag_search, "search_request")
        widgets["contig_select"].param.watch(self.contig_search, "search_request")
        widgets["sample_select"].param.watch(self.sample_search, "search_request")
        widgets["sample_select"].param.watch(self.sample_changed, "value")
        widgets["contig_select"].param.watch(self.contig_changed, "value")
        if widgets["has_mags"]:
            widgets["mag_select"].param.watch(self.mag_changed, "value")
            widgets["contig_select"].param.watch(self.contig_sync_mag, "value")

    def mag_search(self, event: Any) -> None:
        self._search(self.bindings.widgets["mag_select"], self.bindings.compute_mags, event)

    def contig_search(self, event: Any) -> None:
        self._search(self.bindings.widgets["contig_select"], self.bindings.compute_contigs, event)

    def sample_search(self, event: Any) -> None:
        self._search(self.bindings.widgets["sample_select"], self.bindings.compute_samples, event)

    def _search(self, widget: Any, compute: Callable[[str], list[str]], event: Any) -> None:
        if not self.bindings.interaction_lock["locked"]:
            request_nonce, query = decode_search_request(event.new)
            self.bindings.push_completions(widget, compute(query))
            widget.search_result_query = query
            widget.search_result_nonce = request_nonce

    def sample_changed(self, event: Any) -> None:
        if self.bindings.interaction_lock["locked"]:
            return
        with self._held_document():
            valid_samples = self.bindings.compute_samples("")
            self.bindings.push_completions(self.bindings.widgets["sample_select"], valid_samples)
            if event.new and event.new not in valid_samples:
                self.bindings.update_titles()
                return
            if self.bindings.widgets["has_mags"]:
                self.bindings.refresh_mags()
            self.bindings.refresh_contigs()
            self.bindings.update_titles()

    def sync_selected_contig_position(self) -> None:
        bindings = self.bindings
        widgets = bindings.widgets
        if not widgets["has_mags"] or widgets["view_radio"].active != 0:
            return
        mag = widgets["mag_select"].value
        contig = widgets["contig_select"].value
        if not mag or not contig:
            return
        if mag not in widgets["mag_to_contig_offsets"]:
            offset = 0
            offsets = {}
            for member in widgets["mag_to_contigs"].get(mag, []):
                offsets[member] = offset
                offset += widgets["contig_lengths"].get(member, 0)
            widgets["mag_to_contig_offsets"][mag] = offsets
        if contig in widgets["mag_to_contig_offsets"].get(mag, {}):
            offset = widgets["mag_to_contig_offsets"][mag][contig]
            length = widgets["contig_lengths"].get(contig, 0)
            bindings.from_position.value = str(offset + 1)
            bindings.to_position.value = str(offset + length)

    def contig_changed(self, event: Any) -> None:
        if self.bindings.interaction_lock["locked"]:
            return
        with self._held_document():
            widgets = self.bindings.widgets
            if not widgets["has_mags"]:
                self.bindings.refresh_samples()
            self.bindings.update_titles()
            if widgets["has_mags"] and widgets["view_radio"].active == 0:
                mag = widgets["mag_select"].value
                offsets = widgets["mag_to_contig_offsets"].get(mag, {})
                if event.new not in offsets:
                    self.sync_selected_contig_position()
                    offsets = widgets["mag_to_contig_offsets"].get(mag, {})
                offset = offsets.get(event.new)
                length = widgets["contig_lengths"].get(event.new, 0)
                try:
                    start = int(self.bindings.from_position.value)
                    end = int(self.bindings.to_position.value)
                except (TypeError, ValueError):
                    start, end = 0, 0
                contig_start, contig_end = (offset + 1, offset + length) if offset is not None else (1, 0)
                if offset is None or end < contig_start or start > contig_end:
                    self.sync_selected_contig_position()
            elif event.new and event.new in widgets["contig_lengths"]:
                self.bindings.from_position.value = "1"
                self.bindings.to_position.value = str(widgets["contig_lengths"][event.new])
            else:
                self.bindings.from_position.value = "1"
                self.bindings.to_position.value = ""

    def mag_changed(self, event: Any) -> None:
        if self.bindings.interaction_lock["locked"]:
            return
        with self._held_document():
            self.bindings.refresh_contigs()
            self.bindings.refresh_samples()
            self.bindings.update_titles()
            if self.bindings.widgets["view_radio"].active == 0 and event.new:
                total = sum(
                    self.bindings.widgets["contig_lengths"].get(contig, 0)
                    for contig in self.bindings.widgets["mag_to_contigs"].get(event.new, [])
                )
                self.bindings.from_position.value = "1"
                self.bindings.to_position.value = str(total)

    def contig_sync_mag(self, event: Any) -> None:
        if self.bindings.interaction_lock["locked"]:
            return
        with self._held_document():
            parent = self.bindings.widgets["contig_to_mag"].get(event.new) if event.new else None
            if parent and self.bindings.widgets["mag_select"].value != parent:
                self.bindings.widgets["mag_select"].value = parent
            self.bindings.refresh_samples()
            self.bindings.update_titles()

    def subject_scope_changed(self, attr: str, old: int, new: int) -> None:
        with self._held_document():
            self._subject_scope_changed_locked(attr, old, new)

    def _subject_scope_changed_locked(self, attr: str, old: int, new: int) -> None:
        bindings = self.bindings
        widgets = bindings.widgets
        current_order_category = bindings.sample_order_category.value
        if new == 0:
            self.sync_selected_contig_position()
            mag = widgets["mag_select"].value
            contig = widgets["contig_select"].value
            if not (contig and contig in widgets["mag_to_contig_offsets"].get(mag, {})) and mag:
                total = sum(widgets["contig_lengths"].get(c, 0) for c in widgets["mag_to_contigs"].get(mag, []))
                bindings.from_position.value = "1"
                bindings.to_position.value = str(total)
            bindings.sample_current_categories[:] = bindings.sample_mag_categories
            desired_order_category = sample_order_category_for_view(current_order_category, mag_view=True)
        else:
            contig = widgets["contig_select"].value
            if contig and contig in widgets["contig_lengths"]:
                length = widgets["contig_lengths"][contig]
                offset = widgets["mag_to_contig_offsets"].get(widgets["mag_select"].value, {}).get(contig)
                try:
                    start = int(bindings.from_position.value)
                    end = int(bindings.to_position.value)
                except (TypeError, ValueError):
                    start, end = 0, 0
                if offset is None:
                    local_start, local_end = 1, length
                else:
                    local_start, local_end = translate_mag_window_to_contig(
                        start, end, offset=offset, contig_length=length
                    )
                bindings.from_position.value = str(local_start)
                bindings.to_position.value = str(local_end)
            bindings.sample_current_categories[:] = bindings.sample_contig_categories
            desired_order_category = sample_order_category_for_view(current_order_category, mag_view=False)
        bindings.sample_order_category.options = list(bindings.sample_current_categories)
        if desired_order_category not in bindings.sample_current_categories:
            desired_order_category = "Sample" if "Sample" in bindings.sample_current_categories else ""
        if bindings.sample_order_category.value != desired_order_category:
            bindings.sample_order_category.value = desired_order_category
        # A subject-view transition changes the scope of all three selectors.
        # Reproject them together so their option pools and title counts cannot
        # describe different views.
        bindings.refresh_mags()
        bindings.refresh_contigs()
        bindings.refresh_samples()
        bindings.update_titles()

    @contextmanager
    def _held_document(self):
        document = curdoc()
        document.hold("combine")
        self.bindings.interaction_lock["locked"] = True
        try:
            yield document
        finally:
            self.bindings.interaction_lock["locked"] = False
            document.unhold()
