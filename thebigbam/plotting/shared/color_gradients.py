"""Shared numeric gradient interpolation for plot coloring rules."""


def gradient_color(value, value_range, start_color, end_color):
    """Interpolate endpoint colors across a numeric range."""
    low, high = value_range
    fraction = 0.5 if high == low else min(1.0, max(0.0, (value - low) / (high - low)))
    start = tuple(int(start_color[index:index + 2], 16) for index in (1, 3, 5))
    end = tuple(int(end_color[index:index + 2], 16) for index in (1, 3, 5))
    channels = tuple(round(a + (b - a) * fraction) for a, b in zip(start, end))
    return "#" + "".join(f"{channel:02x}" for channel in channels)


def numeric_color_range(values):
    """Return a finite numeric range, ignoring nulls and nonnumeric values."""
    numbers = []
    for value in values:
        try:
            number = float(value)
        except (TypeError, ValueError):
            continue
        if number == number and abs(number) != float("inf"):
            numbers.append(number)
    return (min(numbers), max(numbers)) if numbers else None
