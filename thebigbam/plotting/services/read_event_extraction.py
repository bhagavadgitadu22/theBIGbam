"""Walk CIGAR/MD without expanding matching bases outside the visible window."""

from bisect import bisect_left, bisect_right

from ..models.read_events import ReadEvent
from ..shared.read_cigar import CIGAR
from ..shared.read_md import MD_TOKEN

MAX_EVENT_SEQUENCE = 10_000


def _md_events(md):
    mismatches, deletions, position = {}, {}, 0
    for token in MD_TOKEN.findall(md or ""):
        if token.isdigit():
            position += int(token)
        elif token.startswith("^"):
            deletions[position] = token[1:].upper()
            position += len(token) - 1
        else:
            mismatches[position] = token.upper()
            position += 1
    return mismatches, deletions


def _visible_ranges(position, size, request):
    offset, length = 0, request.contig_length
    low, high = max(0, request.start - 1), min(length, request.end)
    while offset < size:
        local = (position + offset) % length
        width = min(size - offset, length - local)
        a, b = max(local, low), min(local + width, high)
        if a < b:
            yield offset + a - local, offset + b - local
        offset += width


def _point_visible(position, request):
    return (max(0, request.start - 1) <= position % request.contig_length < min(request.end, request.contig_length)
            or (position > 0 and position % request.contig_length == 0 and request.end >= request.contig_length))


def _deletion_bases(deletions, starts, md_position, position, a, b, reference):
    bases = []
    for index in range(a, min(b, a + MAX_EVENT_SEQUENCE)):
        md_index = md_position + index
        token = bisect_right(starts, md_index) - 1
        text_start = starts[token] if token >= 0 else md_index
        text = deletions.get(text_start, "")
        local = md_index - text_start
        base = text[local:local + 1] or reference.base(position + index) or "?"
        bases.append(base)
    return "".join(bases)


def extract_read_events(read, sequence, md, reference, request, budget):
    wanted = request.event_variables
    mismatches, deletions = _md_events(md)
    mismatch_positions = sorted(mismatches)
    deletion_positions = sorted(deletions)
    events, position, query, md_position = [], read.start, 0, 0
    missing, aligned_seen = False, False
    for match in CIGAR.finditer(read.cigar):
        if len(events) >= budget:
            return tuple(events), missing, True
        size, operation = int(match[1]), match[2]
        if operation in "M=X":
            if "mismatches" in wanted and operation != "=":
                unavailable = operation == "M" and md is None and getattr(reference, "unavailable", False)
                missing |= unavailable
                for a, b in (() if unavailable else _visible_ranges(position, size, request)):
                    if operation == "M" and md is not None:
                        indices = [value - md_position for value in mismatch_positions[
                            bisect_left(mismatch_positions, md_position + a):bisect_left(mismatch_positions, md_position + b)
                        ]]
                    else:
                        indices = range(a, b)
                    for index in indices:
                        base = sequence[query + index:query + index + 1].upper()
                        ref = mismatches.get(md_position + index, "") or reference.base(position + index)
                        if operation == "X" or md_position + index in mismatches or (ref and base and base != "=" and base != ref):
                            if len(events) >= budget:
                                return tuple(events), missing, True
                            events.append(ReadEvent("mismatches", position + index, position + index + 1, 1, base, ref))
                            missing |= not bool(base)
                        elif not ref and md is None:
                            missing = True
                            if operation == "M" and getattr(reference, "unavailable", False):
                                break
            position += size
            query += size
            md_position += size
            aligned_seen = True
        elif operation == "I":
            if "insertions" in wanted and _point_visible(position, request):
                events.append(ReadEvent("insertions", position, position, size, sequence[query:query + min(size, MAX_EVENT_SEQUENCE)]))
            query += size
        elif operation == "D":
            if "deletions" in wanted:
                for a, b in _visible_ranges(position, size, request):
                    if len(events) >= budget:
                        return tuple(events), missing, True
                    text = _deletion_bases(deletions, deletion_positions, md_position, position, a, b, reference)
                    missing |= "?" in text
                    events.append(ReadEvent("deletions", position + a, position + b, size, reference=text))
            position += size
            md_position += size
            aligned_seen = True
        elif operation == "N":
            position += size
            aligned_seen = True
        elif operation in "SH":
            kind = "right_clippings" if aligned_seen else "left_clippings"
            if kind in wanted and _point_visible(position, request):
                text = sequence[query:query + min(size, MAX_EVENT_SEQUENCE)] if operation == "S" else ""
                events.append(ReadEvent(kind, position, position, size, text, clipping="soft" if operation == "S" else "hard"))
            if operation == "S":
                query += size
    return tuple(events), missing, False
