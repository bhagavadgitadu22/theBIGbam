"""Prepare annotation dots for the MAG overview track."""

from __future__ import annotations

from ..shared.color_gradients import gradient_color, numeric_color_range


class MagOverviewService:
    def __init__(self, repository, cache=None, profiler=None) -> None:
        self.repository = repository
        self.cache = cache
        self.profiler = profiler

    def _cached(self, key, load):
        if self.cache is None:
            return load()
        hit, value = self.cache.get(key)
        if self.profiler is not None:
            self.profiler.cache("mag_overview", hit)
        if hit:
            return value
        value = load()
        self.cache.put(key, value)
        return value

    def annotation_dots(self, mag_members, color_rules, max_dots=1000):
        if not color_rules or not mag_members:
            return [], [], 0
        # Retrieval depends on membership, not the current visual ordering.
        # Keep the cache reusable when the same MAG is reordered.
        contig_ids = tuple(sorted(cid for cid, _length, _offset in mag_members))
        offsets = {cid: offset for cid, _length, offset in mag_members}
        positions = {
            annotation_id: (contig_id, start, end)
            for annotation_id, contig_id, start, end in self._cached(
                ("mag_overview", "positions", contig_ids),
                lambda: self.repository.annotation_positions(contig_ids),
            )
        }
        colors_by_annotation = {}
        for rule in color_rules:
            try:
                if rule.get("match_mode") == "gradient":
                    key = rule.get("qualifier_key")
                    values = self._cached(
                        ("mag_overview", "values", contig_ids, key),
                        lambda key=key: self.repository.annotation_values(key, contig_ids),
                    )
                    numeric = {}
                    for annotation_id, value in values:
                        value_range = numeric_color_range((value,))
                        if value_range is not None:
                            numeric[annotation_id] = value_range[0]
                    value_range = numeric_color_range(numeric.values())
                    if value_range is None:
                        continue
                    for annotation_id, value in numeric.items():
                        colors_by_annotation.setdefault(
                            annotation_id,
                            gradient_color(
                                value, value_range, rule.get("color", "#2166ac"),
                                rule.get("gradient_end_color", "#b2182b"),
                            ),
                        )
                    continue
                rule_key = (rule.get("qualifier_key"), rule.get("match_mode", "exact"), rule.get("value"))
                matching = self._cached(
                    ("mag_overview", "matches", contig_ids, rule_key),
                    lambda rule=rule: self.repository.matching_annotation_ids(rule, contig_ids),
                )
            except Exception as error:
                print(f"[mag_track] Rule query error (key={rule.get('qualifier_key')}): {error}", flush=True)
                continue
            color = rule.get("color", "#ff0000")
            for annotation_id in matching:
                colors_by_annotation.setdefault(annotation_id, color)
        dots = []
        for annotation_id, color in colors_by_annotation.items():
            if annotation_id not in positions:
                continue
            contig_id, start, end = positions[annotation_id]
            dots.append((offsets.get(contig_id, 0) + (start + end) / 2, color))
        total = len(dots)
        visible = dots if max_dots is None else dots[:max_dots]
        return [x for x, _color in visible], [color for _x, color in visible], total
