"""creTarget CLI banner.

Printed once at startup by `cretarget`. Two rules govern when it appears:

 - It goes to stderr, never stdout, so `cretarget coords-to-bed ... > out.bed`
   stays pipe-safe.
 - By default it appears only when the destination is a terminal. Inside the
   Nextflow pipeline each step runs as its own task with stderr redirected to
   a file, and a genome-wide run has thousands of them -- a banner in every
   `.command.err` is noise, not branding. Pass force=True to override.

Colour follows the NO_COLOR convention (https://no-color.org): the variable
disables colour by its presence, whatever its value.
"""
import os
import sys

_FIG = [
    '                ______                      __ ',
    '  _____________/_  __/___ __________ ____  / /_',
    ' / ___/ ___/ _ \\/ / / __ `/ ___/ __ `/ _ \\/ __/',
    '/ /__/ /  /  __/ / / /_/ / /  / /_/ /  __/ /_  ',
    '\\___/_/   \\___/_/  \\__,_/_/   \\__, /\\___/\\__/  ',
    '                             /____/            ',
]
_ART = [
    '      .-------.    ',
    '     /         v   ',
    ' ---o-------- [##]-',
]
_TAG = "noncoding CRISPRi guide design"


def _isatty(stream) -> bool:
    """isatty() can raise on a detached, closed or exotic stream."""
    try:
        return bool(stream.isatty())
    except Exception:
        return False


def _use_color(stream) -> bool:
    if "NO_COLOR" in os.environ:
        return False
    if os.environ.get("TERM") == "dumb":
        return False
    return _isatty(stream)


def banner(color=None, stream=None) -> str:
    """The banner as a string. `color` defaults to whether `stream` is a
    colour-capable terminal, so a banner written to stdout is not styled
    according to stderr's state."""
    if stream is None:
        stream = sys.stderr
    if color is None:
        color = _use_color(stream)
    red, dim, bold, reset = (
        ("\033[31m", "\033[2m", "\033[1m", "\033[0m") if color else ("", "", "", "")
    )
    art = list(_ART)
    # The red mark is the cut site: the first 'o' on the DNA strand.
    art[2] = art[2].replace("o", red + "o" + reset, 1)
    art[2] += "   " + dim + _TAG + reset
    return bold + "\n".join(_FIG) + reset + "\n\n" + "\n".join(art) + "\n"


def print_banner(file=None, force=False) -> None:
    """Print the banner to `file` (stderr by default).

    Suppressed unless the destination is a terminal, or force=True, or the
    user set CRETARGET_NO_BANNER -- which suppresses it unconditionally.
    """
    if file is None:
        file = sys.stderr
    if os.environ.get("CRETARGET_NO_BANNER"):
        return
    if not force and not _isatty(file):
        return
    try:
        print(banner(stream=file), file=file)
    except Exception:
        # A banner is decoration. It must never be the reason a run fails --
        # a closed stderr or an encoding that cannot represent the box-drawing
        # characters would otherwise raise before any real work starts.
        pass


if __name__ == "__main__":
    print_banner(sys.stdout, force=True)
