#!/usr/bin/env python3
"""
Helper script for invoking clang-format on all the source files in the
project in parallel.
"""

import sys
import subprocess
import xml.etree.ElementTree as ET
import multiprocessing
import argparse
from functools import partial
from concurrent.futures import ThreadPoolExecutor

from util import for_each_source, relative_to_project


def run_clang_format(clang_format, *args):
    """Execute specified clang-format binary with given arguments.

    '-style=file' is appended to the argument list to ensure
    .clang-format file is used. In case of an error (non-zero exit
    status), prints the stderr output and exits with non-zero exit
    status as well.
    """
    call = (clang_format, "-style=file") + args
    try:
        return subprocess.run(call, capture_output=True, text=True, check=True)
    except subprocess.CalledProcessError as exc:
        cmd = " ".join(exc.cmd)
        print(f"{cmd}: non-zero exit status {exc.returncode}")
        print(exc.stderr)
        sys.exit(1)


def check_file(path, tool):
    """Return true if the file is properly formatted.

    Otherwise, print its path (relative to the project directory) to
    stderr and return false.
    """
    res = tool("--output-replacements-xml", path)
    root = ET.fromstring(res.stdout)
    children = list(root)

    is_ok = not children
    if not is_ok:
        relpath = relative_to_project(path)
        print(f"{relpath}", file=sys.stderr)
    return is_ok


def format_file(path, tool):
    """Format file with clang-format and return True."""
    tool("-i", path)
    return True


def run_in_parallel(func, max_threads):
    """Execute function on every source file in parallel.

    Files that are configured by CMake (.in) are omitted.
    """
    with ThreadPoolExecutor(max_threads) as executor:
        future = dict()

        def task(path):
            future[path] = executor.submit(func, path)

        # Do not attempt to format files configured by CMake
        for_each_source(task, excluded=("*.in",))
        results = {path: fut.result() for path, fut in future.items()}

        if not all(results.values()):
            sys.exit(1)


def run(action, args):
    """Execute specified action with configuration given in args."""
    tool = partial(run_clang_format, args.clang_format)
    func = partial(action, tool=tool)
    run_in_parallel(func, max_threads=args.jobs)


def main():
    """Entry point of the script."""
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(title="Actions", dest="action")
    subparsers.add_parser("check", help="check if files are formatted correctly")
    subparsers.add_parser("fix", help="format files")

    max_threads = multiprocessing.cpu_count()
    parser.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=max_threads,
        metavar="N",
        help="number of jobs to run simultaneously",
    )
    parser.add_argument(
        "--clang-format",
        type=str,
        default="clang-format-12",
        metavar="prog",
        help="name or path of clang-format executable to use",
    )
    args = parser.parse_args()

    actions = {
        "check": lambda: run(check_file, args),
        "fix": lambda: run(format_file, args),
        None: parser.print_help,
    }

    action = actions[args.action]
    action()


if __name__ == "__main__":
    main()
