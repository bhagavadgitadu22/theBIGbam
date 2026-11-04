"""Stable alignment identities shared by candidate selection and sampling."""

import hashlib


def read_identity(read):
    return (read.name, read.start, read.end, read.cigar, read.is_read1, read.alignment_class)


def read_score(read):
    payload = repr(read_identity(read)).encode("utf-8")
    return int.from_bytes(hashlib.blake2b(payload, digest_size=8).digest(), "big")
