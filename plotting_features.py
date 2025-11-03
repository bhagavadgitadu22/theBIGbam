from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, HoverTool
import numpy as np
import pysam

# Feature functions
def reduce_position(pos, ref_length):
    return pos % ref_length

def get_coverage(bamfile, reference, ref_length):
    # Preallocate coverage array (integer)
    coverage = np.zeros(ref_length, dtype=np.uint64)

    for read in bamfile.fetch(reference):
        if read.is_unmapped:
            continue

        # For each aligned block (handles indels correctly)
        for start, end in read.get_blocks():
            # if read crosses the circular boundary
            if start < ref_length and end > ref_length:
                coverage[start:ref_length] += 1
                coverage[0:reduce_position(end, ref_length)] += 1
            else:
                coverage[reduce_position(start, ref_length):reduce_position(end, ref_length)] += 1

    return coverage

# Main logic
def smooth_values(coverage, ref_length, window_size):
    # Aggregate by window
    n_windows = (ref_length + window_size - 1) // window_size
    window_cov = np.add.reduceat(coverage, np.arange(0, ref_length, window_size))
    # If genome length not divisible by window size, trim
    window_cov = window_cov[:n_windows]

    # Convert to dict of window -> coverage_depth (sum or mean)
    coverage_dict = {
        i + 1: float(window_cov[i]) / window_size
        for i in range(n_windows)
    }

    return coverage_dict

def make_bokeh_subplot(xx, yy, width, height, color_picked, alpha_picked, title_picked, edge_picked, x_range):
    p = figure(
        width=width,
        height=height,
        title=title_picked,
        x_range=x_range,
        tools="xpan,xwheel_zoom,reset,save"
    )

    # "Step-like" fill (rectangle look)
    if edge_picked:  # We mimic step='pre' by repeating x values to create a step-like area
        step_x = np.repeat(xx, 2)[1:-1]
        step_y = np.repeat(yy, 2)[:-2]
        p.varea(x=step_x, y1=0, y2=step_y, fill_color=color_picked, fill_alpha=alpha_picked)
    else:  # Smooth filled area
        p.varea(x=xx, y1=0, y2=yy, fill_color=color_picked, fill_alpha=alpha_picked)

    # Add hover
    source = ColumnDataSource(data=dict(x=xx, y=yy))
    p.line(x='x', y='y', source=source, line_alpha=0, line_width=1)  # line_alpha=0 makes it invisible
    hover = HoverTool(tooltips=[("Position", "@x"), (title_picked, "@y")], mode='vline')
    p.add_tools(hover)

    # A clean style like your matplotlib setup
    p.yaxis.axis_label = title_picked
    p.yaxis.axis_label_text_font_size = "10pt"
    p.yaxis.axis_label_standoff = 0
    p.ygrid.grid_line_alpha = 0.2
    p.xgrid.visible = False
    p.outline_line_color = None  # hides top/right borders
    p.min_border_left = 40
    p.min_border_right = 10

    return p

def prepare_subplot(bam_file, feature, locus_name, locus_size, max_visible_width, subplot_size, shared_xrange, window_size):
    if feature == "coverage":
        coverage_all_values = get_coverage(bam_file, locus_name, locus_size)
        coverage_averaged = smooth_values(coverage_all_values, locus_size, window_size)
        print(type(coverage_all_values))
        print(type(coverage_averaged))
        print(coverage_all_values[:100])
        print(list(coverage_averaged.values())[:100])

        # Define window midpoints and y-values
        xx = np.arange(window_size/2, locus_size, window_size)
        yy = [coverage_averaged[i + 1] for i in range(len(xx))]

        coverage_fig = make_bokeh_subplot(
            xx=xx,
            yy=yy,
            width=max_visible_width,
            height=subplot_size,
            x_range=shared_xrange,
            color_picked='black',
            alpha_picked=0.5,
            title_picked="Coverage Depth",  # <-- removes the y-axis title
            edge_picked=False
        )
        # Remove the y-axis label
        coverage_fig.yaxis.axis_label = None
        return coverage_fig
    return None

def adding_subplots(mapping_file, features_list, locus_name, locus_size, max_visible_width, subplot_size, shared_xrange, window_size):
    bam_file = pysam.AlignmentFile(mapping_file, "rb")
    bam_references = bam_file.references

    features_list = features_list.split(",")
    subplots = []

    for feature in features_list:
        subplot = prepare_subplot(bam_file, feature, locus_name, locus_size, max_visible_width, subplot_size, shared_xrange, window_size)
        if subplot:
            subplots.append(subplot)

    bam_file.close()
    return subplots