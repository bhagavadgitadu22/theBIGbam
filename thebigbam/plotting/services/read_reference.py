"""Bounded reference-base lookup with circular coordinate support."""

from collections import OrderedDict

import duckdb


class ReadReference:
    def __init__(self, request, cache=None):
        self.source = request.reference_source
        self.contig = request.contig
        self.length = request.contig_length
        self.cache = cache
        self.chunks = OrderedDict()
        self.identity = id(getattr(self.source, "conn", self.source))
        self.unavailable = self.source is None

    def base(self, position):
        if self.unavailable or self.length < 1:
            return ""
        position %= self.length
        start = position // 4096 * 4096
        if start not in self.chunks:
            key = ("read_reference", self.identity, self.contig, start)
            found, sequence = self.cache.get(key) if self.cache is not None else (False, None)
            if not found:
                try:
                    sequence = self.source.sequence(self.contig, start + 1, min(4096, self.length - start)).upper()
                except (duckdb.Error, ValueError, KeyError):
                    sequence = ""
                if self.cache is not None:
                    self.cache.put(key, sequence, retained_bytes=len(sequence) + 128)
            self.chunks[start] = sequence
            if not sequence:
                self.unavailable = True
            while len(self.chunks) > 128:
                self.chunks.popitem(last=False)
        self.chunks.move_to_end(start)
        sequence = self.chunks[start]
        return sequence[position - start:position - start + 1]
