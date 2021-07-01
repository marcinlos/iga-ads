"""
Utility functions for listing, reading and modifying source files.
"""

import os
from pathlib import Path
from fnmatch import fnmatch


# Root source directories
SOURCE_DIRS = ("include", "src", "examples", "tools", "tests")


def project_dir():
    """Return the root directory of the project."""
    script = Path(__file__).resolve()
    return script.parents[1]


def source_roots():
    """
    Return the list of absolute paths of root directories containing
    source and header files.
    """
    root = project_dir()
    return [root / directory for directory in SOURCE_DIRS]


def matches_any(path, rules):
    """Return true if the path matches any of the rules.

    Rules are specified in fnmatch format.
    """
    return any(fnmatch(path, rule) for rule in rules)


def for_each_file(func, directory, suffixes, excluded):
    """Execute function for each matching file in directory."""
    for root, _, files in os.walk(directory):
        for filename in files:
            path = os.path.join(root, filename)
            if filename.endswith(suffixes) and not matches_any(path, excluded):
                func(path)


def for_each_header(func, excluded=()):
    """Execute function for each header file in directory."""
    exts = (".hpp", ".hpp.in")
    for root in source_roots():
        for_each_file(func, root, exts, excluded)


def for_each_source(func, excluded=()):
    """Execute function for each source file (including headers) in directory."""
    exts = (".hpp", ".hpp.in", ".cpp", ".cpp.in")
    for root in source_roots():
        for_each_file(func, root, exts, excluded)


def read_lines(path):
    """Return list of lines comprising specified file.

    Newlines at the end are included.
    """
    with open(path, "r") as file:
        return file.readlines()


def write_lines(path, lines):
    """Write the list of lines to specified file.

    Lines should include newlines.
    """
    with open(path, "w") as file:
        for line in lines:
            file.write(line)
