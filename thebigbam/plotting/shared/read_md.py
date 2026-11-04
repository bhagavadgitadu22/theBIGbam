"""Shared MD tokenization without expanding runs of matching bases."""

import re

MD_TOKEN = re.compile(r"\d+|\^[A-Za-z]+|[A-Za-z]")
