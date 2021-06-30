#!/usr/bin/env python3
"""
Small helper tool for adding, removing and updating the SPDX license
headers to all the C++ source and header files in the project.

Header format:

    SPDX-FileCopyrightText: {from} - {to} {holder <email>}
    SPDX-License-Identifier: {license}
"""

import os
import sys
from pathlib import Path
from fnmatch import fnmatch
import argparse


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


def license_header(year_from, year_to, holder, license_type):
    """Return SPDX license header using specified data."""
    return [
        f"// SPDX-FileCopyrightText: {year_from} - {year_to} {holder}\n",
        f"// SPDX-License-Identifier: {license_type}\n",
    ]


def has_license(lines):
    """Check if first two lines contain a license header."""
    if len(lines) < 2:
        return False

    return (
        "SPDX-FileCopyrightText:" in lines[0] and "SPDX-License-Identifier:" in lines[1]
    )


def check_license_headers():
    """Check if all source files contain valid license headers.

    If files without license header exist, their paths (relative to the
    project directory) are printed and the script exits with non-zero
    exit code.
    """
    invalid = []
    root = project_dir()

    def check(path):
        lines = read_lines(path)
        if not has_license(lines):
            relpath = Path(path).relative_to(root)
            invalid.append(relpath)

    for_each_source(check)

    if invalid:
        for path in invalid:
            print(path)
        sys.exit(1)


def add_header_lines(lines, header):
    """Return list of lines with prepended header."""
    return header + ["\n"] + lines


def remove_header_lines(lines):
    """Return list of lines without the header."""
    lines = lines[2:]
    if lines and lines[0] == "\n":
        lines = lines[1:]
    return lines


def add_header(path, header):
    """Add license header to a file, unless it already has one."""
    lines = read_lines(path)
    if not has_license(lines):
        write_lines(path, add_header_lines(lines, header))


def remove_header(path):
    """Remove license header from a file if it has one."""
    lines = read_lines(path)
    if has_license(lines):
        write_lines(path, remove_header_lines(lines))


def rewrite_header(path, header):
    """Rewrite (or add, if absent) license header to a file."""
    lines = read_lines(path)
    if has_license(lines):
        lines = remove_header_lines(lines)
        lines = add_header_lines(lines, header)
    write_lines(path, lines)


def main():
    """Entry point of the script."""
    year_from = 2015
    year_to = 2021
    email = "marcin.los.91@gmail.com"
    holder = f"Marcin Łoś <{email}>"
    license_type = "MIT"

    header = license_header(year_from, year_to, holder, license_type)

    def with_sources(func):
        return lambda: for_each_source(func)

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(title="Actions", dest="action")
    subparsers.add_parser("check", help="check if license headers are present")
    subparsers.add_parser("add", help="add license headers (if absent)")
    subparsers.add_parser("remove", help="remove license headers")
    subparsers.add_parser("update", help="rewrite license headers")
    args = parser.parse_args()

    actions = {
        "check": check_license_headers,
        "add": with_sources(lambda path: add_header(path, header)),
        "remove": with_sources(remove_header),
        "update": with_sources(lambda path: rewrite_header(path, header)),
        None: parser.print_help
    }

    action = actions[args.action]
    action()


if __name__ == "__main__":
    main()
