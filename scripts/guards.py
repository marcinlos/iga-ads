#!/usr/bin/env python3
"""
Small helper tool for checking and fixing header guards in include
files.

The fixing functionality is very limited - it only works if the guards
are present in exactly the lines they should be present, that is if the
file looks like this:

    // license header spanning
    // first two lines

    #ifndef ...
    #define ...

    <contents>

    #endif ...

In particular, it cannot add missing guards. The final #endif must be
in the last line of the file.

Checking is more flexible and can detect a range of issues.
"""

import sys
import re
import argparse
from util import relative_to_project, read_lines, write_lines, for_each_header


IFNDEF_PATTERN = re.compile(r"#ifndef ([a-zA-Z0-9_]+)\n")
DEFINE_PATTERN = re.compile(r"#define ([a-zA-Z0-9_]+)\n")
ENDIF_PATTERN = re.compile(r"#endif // ([a-zA-Z0-9_]+)\n")


def validate_guards(path, lines):
    """Check if the header guards are present and use correct name.

    The path should be relative to the project root directory.
    """

    class GuardError(Exception):
        """Exception signaling issue with header guards"""

    def ensure(cond, reason):
        if not cond:
            raise GuardError(reason)

    try:
        # 2 lines for license header, 1 empty, 2+1 for header guards
        ensure(len(lines) >= 6, f"File too short ({len(lines)} lines)")

        match_ifndef = IFNDEF_PATTERN.match(lines[3])
        match_define = DEFINE_PATTERN.match(lines[4])
        match_endif = ENDIF_PATTERN.match(lines[-1])
        matches = (match_ifndef, match_define, match_endif)

        ensure(all(matches), f"Missing some lines of guards")

        names = {m.group(1) for m in matches}

        ensure(len(names) == 1, f"Incoherent guard names ({', '.join(names)})")

        name = next(iter(names))
        expected = guard(path)

        ensure(name == expected, f"Invalid guard {name}, should be {expected}")

    except GuardError as ex:
        reason = ex.args[0]
        return reason
    else:
        return None


def check_guards():
    """Check if all header files contain valid header guards.

    If that is not the case, their paths (relative to the project
    directory) are printed together with a message describing the
    issue. The script then exits with non-zero exit code.
    """
    issues = []

    def report(path):
        relpath = relative_to_project(path)
        lines = read_lines(path)

        res = validate_guards(relpath, lines)
        if res is not None:
            issues.append((relpath, res))

    for_each_header(report)

    if issues:
        for path, msg in issues:
            print(f"{path}: {msg}")
        sys.exit(1)


def fix_guards():
    """Rewrite incorrect header guards in all header files."""

    def fix(path):
        relpath = relative_to_project(path)
        lines = read_lines(relpath)

        res = validate_guards(relpath, lines)
        if res is not None:
            macro = guard(relpath)
            lines[3] = f"#ifndef {macro}\n"
            lines[4] = f"#define {macro}\n"
            lines[-1] = f"#endif // {macro}\n"
            write_lines(path, lines)

    for_each_header(fix)


def guard(path):
    """Return header guard macro name for a given path.

    The path should be relative to the project root directory.
    """
    parts = list(path.parts)
    return guard_from_parts(parts[1:])


def guard_from_parts(parts):
    """Combine parts of the path into a guard macro name."""
    name = parts[-1]

    # Remore .in suffix from config files
    if name.endswith(".in"):
        name = name[:-3]

    parts[-1] = name.replace(".", "_")
    return "_".join(part.upper() for part in parts)


def main():
    """Entry point of the script."""
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(title="Actions", dest="action")
    subparsers.add_parser("check", help="check if header guards are correct")
    subparsers.add_parser("fix", help="fix header guards")
    args = parser.parse_args()

    actions = {
        "check": check_guards,
        "fix": fix_guards,
        None: parser.print_help,
    }

    action = actions[args.action]
    action()


if __name__ == "__main__":
    main()
