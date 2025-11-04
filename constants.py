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
def config_feature_subplot(plot_type, color, title, alpha=0.7, size=1):
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
    "reads_starts": config_feature_subplot("bars", "blue", "Reads' Starts"),
    "reads_ends": config_feature_subplot("bars", "blue", "Reads' Ends"),
    "tau": config_feature_subplot("bars", "blue", "Tau"),

    # Misassembly subplots
    "read_lengths": config_feature_subplot("curve", "green", "Read Lengths"),
    "insert_sizes": config_feature_subplot("curve", "green", "Insert Sizes"),
    "bad_orientations": config_feature_subplot("bars", "green", "Bad Orientations"),
    "left_clippings": config_feature_subplot("bars", "purple", "Left Clippings"),
    "right_clippings": config_feature_subplot("bars", "purple", "Right Clippings"),
    "insertions": config_feature_subplot("bars", "red", "Insertions"),
    "deletions": config_feature_subplot("bars", "red", "Deletions"),
    "mismatches": config_feature_subplot("bars", "red", "Mismatches"),
}