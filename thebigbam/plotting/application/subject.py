"""Sample, contig, MAG, search, and genomic-position interaction callbacks."""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from typing import Any, Callable

from bokeh.io import curdoc


CONTIG_TO_MAG_SAMPLE_ORDER_CATEGORY = {
    "Coverage": "MAG coverage",
    "Misassembly": "MAG misassembly",
    "Microdiversity": "MAG microdiversity",
}
MAG_TO_CONTIG_SAMPLE_ORDER_CATEGORY = {
    mag_category: contig_category
    for contig_category, mag_category in CONTIG_TO_MAG_SAMPLE_ORDER_CATEGORY.items()
}


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
        widgets["mag_select"].param.watch(self.mag_search, "search_nonce")
        widgets["contig_select"].param.watch(self.contig_search, "search_nonce")
        widgets["sample_select"].param.watch(self.sample_search, "search_nonce")
        widgets["sample_select"].param.watch(self.sample_changed, "value")
        widgets["contig_select"].param.watch(self.contig_changed, "value")
        if widgets["has_mags"]:
            widgets["mag_select"].param.watch(self.mag_changed, "value")
            widgets["contig_select"].param.watch(self.contig_sync_mag, "value")

    def mag_search(self, event: Any) -> None:
        self._search(self.bindings.widgets["mag_select"], self.bindings.compute_mags)

    def contig_search(self, event: Any) -> None:
        self._search(self.bindings.widgets["contig_select"], self.bindings.compute_contigs)

    def sample_search(self, event: Any) -> None:
        self._search(self.bindings.widgets["sample_select"], self.bindings.compute_samples)

    def _search(self, widget: Any, compute: Callable[[str], list[str]]) -> None:
        if not self.bindings.interaction_lock["locked"]:
            self.bindings.push_completions(widget, compute(widget.search_query))

    def sample_changed(self, event: Any) -> None:
        if self.bindings.interaction_lock["locked"]:
            return
        with self._held_document():
            valid_samples = self.bindings.compute_samples("")
            if event.new and event.new not in valid_samples:
                self.bindings.push_completions(self.bindings.widgets["sample_select"], valid_samples)
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
            desired_order_category = CONTIG_TO_MAG_SAMPLE_ORDER_CATEGORY.get(
                current_order_category, current_order_category
            )
        else:
            contig = widgets["contig_select"].value
            if contig and contig in widgets["contig_lengths"]:
                bindings.from_position.value = "1"
                bindings.to_position.value = str(widgets["contig_lengths"][contig])
            bindings.sample_current_categories[:] = bindings.sample_contig_categories
            desired_order_category = MAG_TO_CONTIG_SAMPLE_ORDER_CATEGORY.get(
                current_order_category, current_order_category
            )
        bindings.sample_order_category.options = list(bindings.sample_current_categories)
        if desired_order_category not in bindings.sample_current_categories:
            desired_order_category = "Sample" if "Sample" in bindings.sample_current_categories else ""
        if bindings.sample_order_category.value != desired_order_category:
            bindings.sample_order_category.value = desired_order_category
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
