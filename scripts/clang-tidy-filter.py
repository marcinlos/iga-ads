#!/usr/bin/env python3
"""Script for deduplicating and filtering clang-tidy output."""

import re
from dataclasses import dataclass
from typing import Tuple
from pathlib import Path
import argparse
from util import relative_to_project


@dataclass(frozen=True)
class WarningInfo:
    """Type used to represent warning occurrence."""

    codes: Tuple[str, ...]  # list of warning names
    file: Path  # path of the file containing the occurrence
    position: Tuple[int, int]  # line and column of occurrence
    text: str  # full text of the warning

    @property
    def line(self) -> int:
        """Line of warning occurrence."""
        return self.position[0]

    @property
    def column(self) -> int:
        """Column of warning occurrence."""
        return self.position[1]

    @property
    def data(self):
        """Complete set of attributes uniquely identifying the occurrence."""
        return (self.codes, self.file, self.position)

    def __eq__(self, other):
        """Equality comparison based on data property."""
        return self.data == other.data

    def __hash__(self):
        """Hash based on data property."""
        return hash(self.data)


# Regex extracting warning names
# e.g. [modernize-avoid-c-arrays,hicpp-move-const-arg]
CODES_PAT = re.compile(r"\[([a-zA-Z,.-]+)\]")

# Regex extracting location information
# e.g. /path/to/ads/include/ads/simulation.hpp:123:31
# If the clang-tidy output is colored, there is one or more ANSI escape
# codes before the actual path, these need to be skipped
LOCATION_PAT = re.compile(r"^(\x1B\[[0-9;]+m)*(.*):(\d+):(\d+):")


def extract_codes(text: str) -> Tuple[str, ...]:
    """Extract names of warnings from full warning text."""
    match = CODES_PAT.search(text)
    if not match:
        raise ValueError("No warning code found")
    return tuple(match.group(1).split(","))


def find_location(text: str) -> re.Match:
    """Find the location information in the full warning text."""
    match = LOCATION_PAT.search(text)
    if not match:
        raise ValueError("No location information found")
    return match


def extract_location(match: re.Match) -> Tuple[Path, Tuple[int, int]]:
    """Create a tuple with warning occurrence location data."""
    path = Path(match.group(2)).resolve()
    path = relative_to_project(path)
    line = int(match.group(3))
    col = int(match.group(4))
    return path, (line, col)


def replace_location(match: re.Match, location: str) -> str:
    """Replace location data in the warning string."""
    text = match.string
    start, end = match.span(2)
    return text[:start] + location + text[end:]


def parse_warning(text: str) -> WarningInfo:
    """Create a WarningInfo object from the full warning text."""
    codes = extract_codes(text)
    match = find_location(text)
    file, pos = extract_location(match)
    new_text = replace_location(match, str(file))
    return WarningInfo(codes, file, pos, new_text)


def starts_new_warning(line) -> bool:
    """Return true if the line starts a new warning."""
    return "warning:" in line


def matches(pattern):
    """Create a predicate that accepts strings matching the pattern."""

    def predicate(text):
        match = re.search(pattern, text)
        return match is not None

    return predicate


def read_warnings(path):
    """Parse clang-tidy output and returns a set of WarningInfo objects."""
    with open(path, "r", encoding="utf-8") as file:
        lines = file.readlines()

    warnings = set()
    current = None

    def end_warning():
        if current:
            warning = parse_warning(current)
            warnings.add(warning)

    def valid(line):
        discard = ("clang-tidy", "Total:")
        return line and not any(s in line for s in discard)

    for line in lines:
        if starts_new_warning(line):
            end_warning()
            current = line
        elif current is not None and valid(line):
            current = current + line

    end_warning()

    return warnings


def make_matcher(args):
    """Return a warning filter based on the CLI arguments."""
    checks = []

    if args.code:
        checks.append(lambda w: any(map(matches(args.code), w.codes)))

    if args.files:
        checks.append(lambda w: matches(args.files)(str(w.file)))

    if args.exclude:
        checks.append(lambda w: not matches(args.exclude)(str(w.file)))

    def satisfied(warning):
        return all(map(lambda f: f(warning), checks))

    return satisfied


def main():
    """Entry point of the script."""
    parser = argparse.ArgumentParser()
    parser.add_argument("file", nargs="?", default="warnings")
    parser.add_argument("--code", "-c")
    parser.add_argument("--files", "-f")
    parser.add_argument("--exclude")
    parser.add_argument("--first-line", action="store_true")
    args = parser.parse_args()

    warnings = read_warnings(args.file)

    pred = make_matcher(args)
    matching = list(filter(pred, warnings))

    def key(warning):
        return (warning.file, warning.line, warning.column)

    matching.sort(key=key)

    for warning in matching:
        if args.first_line:
            text = warning.text.partition("\n")[0] + "\n"
        else:
            text = warning.text
        print(text, end="")

    print(f"\x1B[0mTotal: {len(matching)}")


if __name__ == "__main__":
    main()
