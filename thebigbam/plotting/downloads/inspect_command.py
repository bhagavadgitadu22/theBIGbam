"""Build and present equivalent ``thebigbam inspect`` commands."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Mapping

import panel as pn


@dataclass(frozen=True)
class InspectCommandBindings:
    db_path: str
    widgets: Mapping[str, Any]
    sample_scope: Any
    combined_features: Any | None
    subplot_variables: Mapping[str, list[str]]
    from_position: Any
    to_position: Any
    filtered_samples: Callable[[str], list[str]]
    feature_classifier: Callable[[str], bool] | None = None


class InspectCommandController:
    """Own the download-data button and its command-hint presentation."""

    def __init__(self, bindings: InspectCommandBindings) -> None:
        self.bindings = bindings
        self.pane = pn.pane.HTML("", visible=False, sizing_mode="stretch_width")
        self.button = pn.widgets.Button(name="DOWNLOAD DATA", height=30, button_type="primary", visible=False)
        self.button.on_click(self.show)

    def show(self, event: Any = None) -> None:
        if self.pane.visible:
            self.pane.visible = False
            return
        commands = self.commands()
        if not commands:
            self.pane.object = '<div style="color:#c00;padding:6px">No features or subject selected.</div>'
        else:
            command_html = "".join(
                '<pre style="background:#1a1a2e;color:#e0e0e0;padding:8px 12px;border-radius:4px;'
                "font-size:13px;white-space:pre-wrap;word-break:break-all;margin:4px 0;"
                'user-select:all;cursor:pointer" title="Click to select, then Ctrl+C to copy">'
                f"{command}</pre>"
                for command in commands
            )
            self.pane.object = (
                '<div style="padding:6px 0"><div style="font-size:12px;color:#888;margin-bottom:4px">'
                "Run in your terminal to download data as TSV (click command to select):</div>"
                f"{command_html}</div>"
            )
        self.pane.visible = True

    def commands(self) -> list[str]:
        bindings = self.bindings
        widgets = bindings.widgets
        has_samples = widgets["has_samples"]
        all_samples = bindings.sample_scope.active == 1 if has_samples else False
        mag_view = widgets["has_mags"] and widgets["view_radio"].active == 0
        features = self._selected_features(all_samples)
        if not features:
            return []

        if mag_view:
            subject = widgets["mag_select"].value
            if not subject:
                return []
            length = sum(widgets["contig_lengths"].get(c, 0) for c in widgets["mag_to_contigs"].get(subject, []))
            subject_arg = f"--mag {subject}"
            samples = [widgets["sample_select"].value] if has_samples and widgets["sample_select"].value else []
        else:
            subject = widgets["contig_select"].value
            if not subject:
                return []
            length = widgets["contig_lengths"].get(subject, 0)
            subject_arg = f"--contig {subject}"
            if not has_samples:
                samples = []
            elif all_samples:
                samples = bindings.filtered_samples(subject)
            else:
                selected = widgets["sample_select"].value
                samples = [selected] if selected else []

        region = self._region_argument(length)
        classifier = bindings.feature_classifier
        if classifier is None:
            from thebigbam.database.blob_decoder import is_contig_blob_feature

            classifier = is_contig_blob_feature
        contig_features = [feature for feature in features if classifier(feature)]
        sample_features = [feature for feature in features if not classifier(feature)]
        prefix = f"thebigbam inspect -d {bindings.db_path} {subject_arg}"
        commands = []
        if contig_features:
            commands.append(f"{prefix} --feature {','.join(contig_features)}{region} > output.tsv")
        if sample_features:
            sample_arg = f" --sample {','.join(samples)}" if samples else " --sample &lt;SAMPLE&gt;"
            commands.append(f"{prefix}{sample_arg} --feature {','.join(sample_features)}{region} > output.tsv")
        return commands

    def _selected_features(self, all_samples: bool) -> list[str]:
        bindings = self.bindings
        features = []
        combined = bindings.combined_features
        if combined is not None:
            for index in combined.active:
                label = combined.labels[index]
                features.extend(bindings.subplot_variables.get(label, [label]))
        groups = bindings.widgets["variables_widgets_all" if all_samples else "variables_widgets_one"]
        for group in groups:
            active = group.active[-1:] if all_samples else group.active
            for index in active:
                label = group.labels[index]
                features.extend(bindings.subplot_variables.get(label, [label]))
            if all_samples and active:
                break
        return features

    def _region_argument(self, length: int) -> str:
        try:
            start = int(self.bindings.from_position.value) if self.bindings.from_position.value.strip() else 1
            end = int(self.bindings.to_position.value) if self.bindings.to_position.value.strip() else length
        except ValueError:
            start, end = 1, length
        if length > 0:
            start = max(1, min(start, length))
            end = max(1, min(end, length))
        return f" --region {start}-{end}" if length > 0 and start < end and (start > 1 or end < length) else ""
