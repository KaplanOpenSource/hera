"""The filter that pins a chosen list of packages from the full requirements."""
import os
import re
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from minimal_requirements import minimal_requirements, package_name

FULL = """\
# Core dependencies
fastapi>=0.115.0
uvicorn[standard]>=0.30.0
Deprecated==1.2.18
dask[dataframe]==2024.8.0

Unum==4.2.1
matplotlib==3.9.4
matplotlib-inline==0.1.7
"""


class TestPackageName:
    def test_a_plain_name(self):
        assert package_name("fastapi\n") == "fastapi"

    def test_drops_the_version(self):
        assert package_name("Unum==4.2.1") == "unum"

    def test_drops_the_extras(self):
        assert package_name("uvicorn[standard]>=0.30.0") == "uvicorn"

    def test_dots_and_underscores_become_dashes(self):
        assert package_name("matplotlib_inline==0.1.7") == "matplotlib-inline"
        assert package_name("ruamel.yaml==0.18.6") == "ruamel-yaml"

    def test_drops_an_environment_marker(self):
        assert package_name("tomli==2.0.1; python_version < '3.11'") == "tomli"

    def test_drops_a_trailing_comment(self):
        assert package_name("numpy==1.26.4  # pinned") == "numpy"

    def test_a_blank_line_names_nothing(self):
        assert package_name("   \n") == ""

    def test_a_comment_keeps_its_hash(self):
        assert package_name("# Core dependencies") == "#"

    def test_an_option_keeps_its_dash(self):
        assert package_name("-r other.txt") == "-r"


@pytest.fixture
def full(tmp_path):
    path = tmp_path / "requirements.txt"
    path.write_text(FULL)
    return str(path)


def chosen(tmp_path, text):
    """A chosen list holding the given lines."""
    path = tmp_path / "chosen.txt"
    path.write_text(text)
    return str(path)


def run(full, tmp_path, text, capsys):
    """The lines the filter prints for the given chosen list."""
    minimal_requirements(full, chosen(tmp_path, text))
    return capsys.readouterr().out.splitlines()


def test_takes_the_pin_of_a_name(full, tmp_path, capsys):
    assert run(full, tmp_path, "Unum\n", capsys) == ["Unum==4.2.1"]


def test_keeps_the_extras_the_full_file_has(full, tmp_path, capsys):
    assert run(full, tmp_path, "dask[dataframe]\n", capsys) == ["dask[dataframe]==2024.8.0"]


def test_a_name_without_extras_finds_the_line_with_them(full, tmp_path, capsys):
    assert run(full, tmp_path, "uvicorn\n", capsys) == ["uvicorn[standard]>=0.30.0"]


def test_the_case_does_not_matter(full, tmp_path, capsys):
    lines = run(full, tmp_path, "UNUM\ndeprecated\n", capsys)

    assert lines == ["Unum==4.2.1", "Deprecated==1.2.18"]


def test_a_dash_is_the_same_as_a_dot_or_an_underscore(full, tmp_path, capsys):
    assert run(full, tmp_path, "matplotlib_inline\n", capsys) == ["matplotlib-inline==0.1.7"]


def test_a_shorter_name_does_not_take_a_longer_one(full, tmp_path, capsys):
    assert run(full, tmp_path, "matplotlib\n", capsys) == ["matplotlib==3.9.4"]


def test_a_greater_than_pin_comes_as_it_is(full, tmp_path, capsys):
    assert run(full, tmp_path, "fastapi\n", capsys) == ["fastapi>=0.115.0"]


def test_comments_blanks_and_options_are_skipped(full, tmp_path, capsys):
    assert run(full, tmp_path, "# a comment\n\n-r other.txt\nUnum\n", capsys) == ["Unum==4.2.1"]


def test_the_order_is_the_chosen_order(full, tmp_path, capsys):
    lines = run(full, tmp_path, "Unum\nfastapi\n", capsys)

    assert lines == ["Unum==4.2.1", "fastapi>=0.115.0"]


def test_a_name_the_full_file_lacks_comes_as_it_is(full, tmp_path, capsys):
    assert run(full, tmp_path, "argcomplete>=3.1.0\n", capsys) == ["argcomplete>=3.1.0"]


def test_a_bare_name_the_full_file_lacks_comes_bare(full, tmp_path, capsys):
    assert run(full, tmp_path, "argcomplete\n", capsys) == ["argcomplete"]


def test_the_full_file_gives_every_pin_in_it(full, tmp_path, capsys):
    minimal_requirements(full, full)

    assert capsys.readouterr().out.splitlines() == [
        "fastapi>=0.115.0",
        "uvicorn[standard]>=0.30.0",
        "Deprecated==1.2.18",
        "dask[dataframe]==2024.8.0",
        "Unum==4.2.1",
        "matplotlib==3.9.4",
        "matplotlib-inline==0.1.7",
    ]


REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

class TestTheRepoList:
    """The chosen list must find all its pins, or a build installs whatever pip likes."""

    def test_every_name_is_pinned_in_requirements(self, capsys):
        full = os.path.join(REPO, "requirements.txt")

        minimal_requirements(full, os.path.join(REPO, "requirements-minimal.txt"))

        unpinned = [line for line in capsys.readouterr().out.splitlines()
                    if not re.search(r"[=<>~]", line)]
        assert unpinned == []
