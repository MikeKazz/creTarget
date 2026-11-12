"""The banner is decoration, so these tests are about what it must never do:
contaminate stdout, appear in non-interactive logs, or raise.
"""
import io
import sys

import pytest

from cretarget import banner as b


class _Tty(io.StringIO):
    def isatty(self):
        return True


class _NotTty(io.StringIO):
    def isatty(self):
        return False


class _Hostile(io.StringIO):
    """isatty() raising is not hypothetical: it happens on streams detached
    by a job scheduler."""

    def isatty(self):
        raise OSError("detached")


def test_art_is_rectangular():
    assert len({len(line) for line in b._FIG}) == 1
    assert len({len(line) for line in b._ART}) == 1


def test_no_color_when_not_a_terminal():
    out = b.banner(stream=_NotTty())
    assert "\033" not in out


def test_color_when_a_terminal(monkeypatch):
    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.delenv("TERM", raising=False)
    assert "\033[31m" in b.banner(stream=_Tty())


@pytest.mark.parametrize("value", ["", "1", "anything"])
def test_no_color_env_disables_color_whatever_its_value(monkeypatch, value):
    """no-color.org specifies presence, not truthiness."""
    monkeypatch.setenv("NO_COLOR", value)
    assert "\033" not in b.banner(stream=_Tty())


def test_dumb_terminal_gets_no_color(monkeypatch):
    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.setenv("TERM", "dumb")
    assert "\033" not in b.banner(stream=_Tty())


def test_suppressed_when_destination_is_not_a_terminal(monkeypatch):
    """The reason this matters: a genome-wide run has thousands of tasks, each
    with stderr redirected to a file."""
    monkeypatch.delenv("CRETARGET_NO_BANNER", raising=False)
    f = _NotTty()
    b.print_banner(f)
    assert f.getvalue() == ""


def test_force_overrides_the_terminal_check(monkeypatch):
    monkeypatch.delenv("CRETARGET_NO_BANNER", raising=False)
    f = _NotTty()
    b.print_banner(f, force=True)
    assert "creTarget" not in f.getvalue()  # it is ASCII art, not literal text
    assert "noncoding CRISPRi guide design" in f.getvalue()


def test_env_var_suppresses_even_with_force(monkeypatch):
    monkeypatch.setenv("CRETARGET_NO_BANNER", "1")
    f = _Tty()
    b.print_banner(f, force=True)
    assert f.getvalue() == ""


def test_never_raises_on_a_hostile_stream(monkeypatch):
    monkeypatch.delenv("CRETARGET_NO_BANNER", raising=False)
    b.print_banner(_Hostile())          # must not raise
    b.print_banner(_Hostile(), force=True)


def test_cli_leaves_stdout_clean(monkeypatch, capsys, tmp_path):
    """`cretarget coords-to-bed ... > out.bed` must stay pipe-safe."""
    from cretarget import cli

    monkeypatch.delenv("CRETARGET_NO_BANNER", raising=False)
    src = tmp_path / "coords.tsv"
    src.write_text("region_id\nchr1:1000-2000\n")
    out = tmp_path / "out.bed"

    cli.main(["coords-to-bed", "--input", str(src),
              "--column", "region_id", "--out", str(out)])

    captured = capsys.readouterr()
    assert "\033" not in captured.out
    assert "/____/" not in captured.out
    assert out.read_text().startswith("chr1\t")
