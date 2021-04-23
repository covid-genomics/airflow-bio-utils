from __future__ import with_statement

import mmap
from typing import Optional


def get_file_lines_count(filename: str, text: Optional[str] = None) -> int:
    """
    Count lines in file.
    This is the fastest implementation possible in Python using mmap.

    Optionally function takes a text parameter.
    If that parameter is not None, then only lines containing this text are counted.

    :param filename: Path to the input file
    :param text: Optional text to match in counted lines
    :return: Number of counted lines
    """
    f = open(filename, "r+")
    buf = mmap.mmap(f.fileno(), 0)
    lines = 0
    readline = buf.readline
    if text is None:
        while readline():
            lines += 1
    else:
        text = text.encode()
        for line in buf:
            if text in line:
                lines += 1
    return lines
