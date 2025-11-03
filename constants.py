# Define function-to-color mapping
# Use the color scheme from pharokka
PHAROKKA_COLORS = {
    "vfdb_card": "#FF0000",
    "unknown function": "#AAAAAA",
    "other": "#4deeea",
    "tail": "#74ee15",
    "transcription regulation": "#ffe700",
    "dna, rna and nucleotide metabolism": "#f000ff",
    "lysis": "#001eff",
    "moron, auxiliary metabolic gene and host takeover": "#8900ff",
    "integration and excision": "#E0B0FF",
    "head and packaging": "#ff008d",
    "connector": "#5A5A5A",
}

# Define subplots characteristics
def config_feature_subplot(plot_type, color, title, alpha=0.8, size=2):
    if plot_type == "bars":
        alpha = 0.5
        size = 3
    return {
        "type_picked": plot_type,
        "color_picked": color,
        "alpha_picked": alpha,
        "size_picked": size,
        "title_picked": title
    }

FEATURE_SUBPLOTS = {
    "coverage": config_feature_subplot("curve", "black", "Coverage Depth"),

    # Starts subplots
    "coverage_reduced": config_feature_subplot("curve", "black", "Coverage Depth (only reads starting and ending with a match)"),
    "reads_starts": config_feature_subplot("bars", "blue", "Reads' Starts (more than deviation_factor*std away from the mean)"),
    "reads_ends": config_feature_subplot("bars", "blue", "Reads' Ends (more than deviation_factor*std away from the mean)"),
    "tau": config_feature_subplot("bars", "blue", "Tau (more than deviation_factor*std away from the mean)"),

    # Misassembly subplots
    "read_lengths": config_feature_subplot("bars", "green", "Read Lengths (more than deviation_factor*std away from the mean)"),
    "insert_sizes": config_feature_subplot("bars", "green", "Insert Sizes (more than deviation_factor*std away from the mean)"),
    "bad_orientations": config_feature_subplot("bars", "green", "Bad Orientations (more than deviation_factor*std away from the mean)"),
    "left_clippings": config_feature_subplot("bars", "purple", "Left Clippings (more than deviation_factor*std away from the mean)"),
    "right_clippings": config_feature_subplot("bars", "purple", "Right Clippings (more than deviation_factor*std away from the mean)"),
    "insertions": config_feature_subplot("bars", "red", "Insertions (more than deviation_factor*std away from the mean)"),
    "deletions": config_feature_subplot("bars", "red", "Deletions (more than deviation_factor*std away from the mean)"),
    "mismatches": config_feature_subplot("bars", "pink", "Mismatches (more than deviation_factor*std away from the mean)"),
}