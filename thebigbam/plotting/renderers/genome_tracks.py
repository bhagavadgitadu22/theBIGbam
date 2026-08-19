"""Bokeh-only renderers for annotation and sequence track data."""

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from bokeh.models import (
    BoxZoomTool,
    ColumnDataSource,
    HoverTool,
    NumeralTickFormatter,
    Range1d,
    SaveTool,
    TapTool,
    WheelZoomTool,
)
from bokeh.plotting import figure
from dna_features_viewer import BiopythonTranslator

from ..services.genome_tracks import GeneMapData, NucleotideSequenceData, TranslationData


class _Translator(BiopythonTranslator):
    def compute_feature_color(self, feature):
        return feature.qualifiers.get("_custom_color") or ("#cccccc" if feature.type == "CDS" else "#000000")

    def compute_feature_label(self, feature):
        return None

    def compute_feature_html(self, feature):
        key = feature.qualifiers.get("_tooltip_key", "product")
        return str(feature.qualifiers.get(key) or feature.type)


def _finish_track(plot):
    plot.toolbar.logo = None
    plot.xaxis.formatter = NumeralTickFormatter(format="0,0")
    wheel = WheelZoomTool(dimensions="width")
    plot.add_tools(wheel, BoxZoomTool(dimensions="width"))
    plot.toolbar.active_scroll = wheel
    return plot


class GeneMapRenderer:
    def render(self, data: GeneMapData, height, x_range, figure_width=30):
        features = []
        for item in data.features:
            qualifiers = {key: value for key, value in item.qualifiers.items() if value is not None}
            if item.color:
                qualifiers["_custom_color"] = item.color
            if item.tooltip_key:
                qualifiers["_tooltip_key"] = item.tooltip_key
            try:
                location = FeatureLocation(item.start - 1, item.end, strand=item.strand)
            except ValueError:
                continue
            features.append(SeqFeature(location=location, type=item.feature_type, qualifiers=qualifiers))
        if not features:
            return None
        record = SeqRecord(Seq("N" * data.subject_length), id=data.subject_name, features=features)
        plot = _Translator().translate_record(record).plot_with_bokeh(figure_width=figure_width, figure_height=40)
        plot.height = height
        plot.tools = [tool for tool in plot.tools if not isinstance(tool, TapTool)]
        plot.add_tools(SaveTool())
        plot.x_range = x_range
        return _finish_track(plot)


class NucleotideRenderer:
    COLORS = {"A": "#d62728", "T": "#2ca02c", "G": "#ff7f0e", "C": "#1f77b4"}

    def render(self, data: NucleotideSequenceData, height, x_range):
        if not data.positions:
            return None
        positions = list(data.positions)
        letters = list(data.nucleotides)
        source = ColumnDataSource(
            data={
                "left": [p - 0.5 for p in positions],
                "right": [p + 0.5 for p in positions],
                "bottom": [0] * len(positions),
                "top": [1] * len(positions),
                "color": [self.COLORS.get(nt.upper(), "#999999") for nt in letters],
                "nucleotide": letters,
                "position": positions,
                "x_center": [float(p) for p in positions],
                "y_center": [0.5] * len(positions),
            }
        )
        plot = figure(height=height, x_range=x_range, y_range=Range1d(0, 1), tools="xpan,reset,save")
        plot.quad(left="left", right="right", bottom="bottom", top="top", color="color", source=source, line_color=None)
        plot.text(
            x="x_center",
            y="y_center",
            text="nucleotide",
            source=source,
            text_align="center",
            text_baseline="middle",
            text_font_size="8pt",
            text_color="white",
            text_font_style="bold",
        )
        plot.add_tools(HoverTool(tooltips=[("Position", "@position{0,0}"), ("Nucleotide", "@nucleotide")]))
        plot.yaxis.visible = False
        plot.grid.visible = False
        plot.outline_line_color = None
        return _finish_track(plot)


def _codon_quads(segments, strand, index):
    ordered = segments if strand >= 0 else tuple(reversed(segments))
    wanted_start, wanted_end, cursor = index * 3, index * 3 + 3, 0
    result = []
    for start, end in ordered:
        segment_end = cursor + end - start + 1
        lo, hi = max(wanted_start, cursor), min(wanted_end, segment_end)
        if lo < hi:
            if strand >= 0:
                left, right = start + lo - cursor, start + hi - cursor - 1
            else:
                left, right = end - (hi - cursor - 1), end - (lo - cursor)
            result.append((left - 0.5, right + 0.5))
        cursor = segment_end
        if cursor >= wanted_end:
            break
    return result


class TranslationRenderer:
    def render(self, data: TranslationData, start, end, height, x_range):
        if not data.entries:
            return None
        lane_ends, lanes = [], {}
        for index, entry in sorted(enumerate(data.entries), key=lambda pair: pair[1].start):
            lane = next((i for i, lane_end in enumerate(lane_ends) if lane_end < entry.start), None)
            if lane is None:
                lane = len(lane_ends)
                lane_ends.append(entry.end)
            else:
                lane_ends[lane] = entry.end
            lanes[index] = lane
        columns = {
            key: []
            for key in (
                "left",
                "right",
                "bottom",
                "top",
                "color",
                "amino_acid",
                "amino_acid_name",
                "aa_label",
                "codon",
                "position_start",
                "position_end",
                "x_center",
                "y_center",
            )
        }
        count = len(lane_ends)
        for index, entry in enumerate(data.entries):
            for aa_index, amino_acid in enumerate(entry.protein_sequence):
                codon = entry.nucleotide_sequence[aa_index * 3 : aa_index * 3 + 3].upper()
                info = data.codons.get(codon, (amino_acid, "Unknown", "???", "#999999"))
                for left, right in _codon_quads(entry.segments, entry.strand, aa_index):
                    if right < start or left > end:
                        continue
                    top, bottom = 1 - lanes[index] / count, 1 - (lanes[index] + 1) / count
                    values = (
                        left,
                        right,
                        bottom,
                        top,
                        info[3],
                        amino_acid,
                        info[1],
                        info[2],
                        codon,
                        int(left + 0.5),
                        int(right - 0.5),
                        (left + right) / 2,
                        (bottom + top) / 2,
                    )
                    for key, value in zip(columns, values):
                        columns[key].append(value)
        if not columns["left"]:
            return None
        source = ColumnDataSource(data=columns)
        plot = figure(height=height, x_range=x_range, y_range=Range1d(0, 1), tools="xpan,reset,save")
        plot.quad(left="left", right="right", bottom="bottom", top="top", color="color", source=source, line_color=None)
        plot.text(
            x="x_center",
            y="y_center",
            text="aa_label",
            source=source,
            text_align="center",
            text_baseline="middle",
            text_font_size="7pt",
            text_color="white",
            text_font_style="bold",
        )
        plot.add_tools(
            HoverTool(
                tooltips=[
                    ("Position", "@position_start{0,0}-@position_end{0,0}"),
                    ("Codon", "@codon"),
                    ("Amino acid", "@amino_acid (@amino_acid_name)"),
                ]
            )
        )
        plot.yaxis.visible = False
        plot.grid.visible = False
        plot.outline_line_color = None
        return _finish_track(plot)
