"""Count exact matches from CIGAR and standard edit-distance/MD tags."""

from ..shared.read_md import MD_TOKEN


def exact_matching_bases(read):
    counts = read.get_cigar_stats()[0]
    m, insertions, deletions, soft, hard, equal, different = (int(counts[i]) for i in (0, 1, 2, 4, 5, 7, 8))
    read_length = m + insertions + soft + hard + equal + different
    if not m:
        return equal, read_length
    aligned = m + equal + different
    if read.has_tag("NM"):
        try:
            edits = read.get_tag("NM")
            mismatch_count = int(edits) - insertions - deletions
            if float(edits).is_integer() and different <= mismatch_count <= m + different:
                return aligned - mismatch_count, read_length
        except (TypeError, ValueError, OverflowError):
            pass
    if read.has_tag("MD"):
        mismatch_count = reference_length = deleted = 0
        for match in MD_TOKEN.finditer(str(read.get_tag("MD"))):
            token = match[0]
            if token.isdigit():
                reference_length += int(token)
            elif token.startswith("^"):
                deleted += len(token) - 1
                reference_length += len(token) - 1
            else:
                mismatch_count += 1
                reference_length += 1
        if reference_length == aligned + deletions and deleted == deletions and different <= mismatch_count <= m + different:
            return aligned - mismatch_count, read_length
    return None, read_length
