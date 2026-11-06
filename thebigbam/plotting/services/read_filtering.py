"""Compile shared AND/OR expressions into alignment predicates before sampling."""

import math
import operator
from time import perf_counter

from ..models.read_filters import (
    BOOLEAN_READ_FILTERS,
    READ_FILTER_CHOICES,
    READ_FILTER_FIELDS,
    READ_FILTER_LABELS,
    read_filter_operators,
)

COMPARISONS = {">": operator.gt, "<": operator.lt, "=": operator.eq, "!=": operator.ne}


def _boolean(value):
    if isinstance(value, bool):
        return value
    if isinstance(value, str) and value.upper() in {"TRUE", "FALSE"}:
        return value.upper() == "TRUE"
    raise ValueError("Boolean read filters require TRUE or FALSE.")


def validate_read_filters(expression):
    if any(value not in {"AND", "OR"} for value in expression.connectors):
        raise ValueError("Read filter blocks must use AND or OR.")
    for section in expression.sections:
        if any(value not in {"AND", "OR"} for value in section.connectors):
            raise ValueError("Read filter rows must use AND or OR.")
        for rule in section.predicates:
            label = READ_FILTER_LABELS.get(rule.column, rule.column)
            if rule.category != "Reads" or rule.column not in READ_FILTER_FIELDS:
                raise ValueError(f"Unknown read filter: {rule.column}.")
            if rule.operator not in read_filter_operators(rule.column):
                raise ValueError(f"Unsupported operator for {label}: {rule.operator}.")
            if rule.column in BOOLEAN_READ_FILTERS:
                _boolean(rule.value)
                continue
            if rule.column in READ_FILTER_CHOICES:
                if rule.value not in READ_FILTER_CHOICES[rule.column]:
                    choices = " or ".join(READ_FILTER_CHOICES[rule.column])
                    raise ValueError(f"{label} must be {choices}.")
                continue
            try:
                number = float(rule.value)
            except (TypeError, ValueError) as error:
                raise ValueError(f"{label} requires a number.") from error
            if not math.isfinite(number):
                raise ValueError(f"{label} requires a finite number.")
            if rule.column == "Read % aligned" and not 0 <= number <= 100:
                raise ValueError("Read % aligned must be between 0 and 100.")
            if rule.column == "Alignment includes" and (number < 1 or not number.is_integer()):
                raise ValueError("Alignment includes requires a positive integer position.")


def _has_mate(read):
    return bool(read.is_paired and not read.flag & 8 and read.mate_start is not None and read.mate_start >= 0)


def _non_inward(read, length, circular):
    if not _has_mate(read) or not read.mate_same_reference:
        return False
    reverse, mate_reverse = bool(read.flag & 16), bool(read.flag & 32)
    if read.start == read.mate_start:
        return reverse == mate_reverse
    read_is_left = (read.mate_start - read.start) % length <= length / 2 if circular else read.start < read.mate_start
    inward = (not reverse and mate_reverse) if read_is_left else (reverse and not mate_reverse)
    return not inward


def _value(read, field, length, circular):
    if field == "Read length":
        return read.lengths[0]
    if field == "Alignment length":
        return read.lengths[1]
    if field == "Read % aligned":
        return read.alignment_percentage
    if field == "Alignment start":
        return read.start + 1
    if field == "Alignment end":
        return read.end
    if field == "MAPQ":
        return read.mapq
    if field == "Strand":
        return "Forward" if read.forward else "Reverse"
    if field == "Alignment class":
        return read.alignment_class.title()
    if field == "Insert sizes":
        return abs(read.template_length) if _has_mate(read) and read.mate_same_reference and read.template_length else None
    if field == "Has mate mapped":
        return _has_mate(read)
    if field == "Non-inward pair":
        return _non_inward(read, length, circular)
    return bool(circular and read.end > length)


def _includes(read, position, length, circular):
    if not 1 <= position <= length:
        return False
    if not circular:
        return read.start < position <= read.end
    if read.end - read.start >= length:
        return True
    start, end = read.start % length + 1, (read.end - 1) % length + 1
    return start <= position <= end if start <= end else position >= start or position <= end


def _fold(functions, connectors, read):
    result = functions[0](read)
    for connector, function in zip(connectors, functions[1:]):
        result = (result or function(read)) if connector == "OR" else (result and function(read))
    return result


def compile_read_filters(expression, length, circular):
    validate_read_filters(expression)
    sections, connectors = [], []
    for index, section in enumerate(expression.sections):
        predicates = []
        for rule in section.predicates:
            if rule.column == "Alignment includes":
                position, negate = int(float(rule.value)), rule.operator == "!="
                predicates.append(lambda read, position=position, negate=negate: _includes(read, position, length, circular) != negate)
            else:
                expected = (
                    _boolean(rule.value)
                    if rule.column in BOOLEAN_READ_FILTERS
                    else rule.value if rule.column in READ_FILTER_CHOICES
                    else float(rule.value)
                )
                compare = COMPARISONS[rule.operator]

                def predicate(read, field=rule.column, expected=expected, compare=compare):
                    actual = _value(read, field, length, circular)
                    return actual is not None and compare(actual, expected)

                predicates.append(predicate)
        if predicates:
            if sections:
                connectors.append(expression.connectors[index - 1])
            sections.append(lambda read, predicates=predicates, joins=section.connectors: _fold(predicates, joins, read))
    return (lambda read: _fold(sections, connectors, read)) if sections else (lambda read: True)


def filtered_candidates(candidates, expression, length, circular, profiler=None):
    matches = compile_read_filters(expression, length, circular)
    for read in candidates:
        before = perf_counter() if profiler is not None else 0
        accepted = matches(read)
        if profiler is not None:
            profiler.add_time("read_filtering", perf_counter() - before)
            profiler.count("read_filter_examined")
            profiler.count("read_filter_matched", int(accepted))
        if accepted:
            yield read
