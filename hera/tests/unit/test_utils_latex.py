"""bibItem / bibtexFile: converting .bbl bibliographies to Hebrew LaTeX.

The conversion rules are stated in bibItem.convert's own comments:

  * wrap English words in \\L{}
  * replace "and" with a leading 'ו' on the next word (dropping the comma
    that preceded it)
  * wrap numbers in $...$
  * leave LaTeX commands alone

Assertions below follow those rules.
"""
import io

import pytest

from hera.utils.latex import bibItem, bibtexFile


@pytest.mark.unit
class TestHebrewDetection:
    def test_detects_hebrew_in_any_line(self):
        assert bibItem(["English only.", "יוסף כוהן."]).isHebrew is True

    def test_pure_english_is_not_hebrew(self):
        assert bibItem(["John Smith.", "Wind speed."]).isHebrew is False

    def test_empty_item_is_not_hebrew(self):
        assert bibItem([]).isHebrew is False

    def test_line_level_detection(self):
        item = bibItem([])
        assert item.isHebrewText("שלום") is True
        assert item.isHebrewText("hello") is False
        assert item.isHebrewText("123 !@#") is False

    def test_unicode_name_of_a_nameless_character(self):
        """Control characters have no Unicode name; 'NONE' must come back."""
        item = bibItem([])
        assert item._getUnicodeName("A") == "LATIN CAPITAL LETTER A"
        assert item._getUnicodeName("\x00") == "NONE"


@pytest.mark.unit
class TestNonHebrewConversion:
    def test_prefixed_with_unsethebrew_and_left_verbatim(self):
        converted = bibItem(["John Smith. Wind speed."]).convert()
        assert converted == "\\unsethebrew\nJohn Smith. Wind speed."

    def test_no_latin_wrapping_is_applied(self):
        converted = bibItem(["John Smith."]).convert()
        assert "\\L{" not in converted


@pytest.mark.unit
class TestHebrewConversion:
    def test_prefixed_with_sethebrew(self):
        assert bibItem(["יוסף כוהן."]).convert().startswith("\\sethebrew\n")

    def test_latin_words_are_wrapped(self):
        converted = bibItem(["יוסף Cohen כתב"]).convert()
        assert "\\L{Cohen}" in converted

    def test_numbers_are_wrapped_in_math_mode(self):
        converted = bibItem(["מהירות 42 קמ"]).convert()
        assert "$42$" in converted

    def test_and_becomes_a_prefixed_vav(self):
        """"and" joins the next word as 'ו', and the preceding comma goes."""
        converted = bibItem(["א, and ב"]).convert()
        assert converted == "\\sethebrew\nא וב"

    def test_latex_commands_are_untouched(self):
        converted = bibItem(["\\newblock שלום"]).convert()
        assert "\\newblock" in converted
        assert "\\L{newblock}" not in converted

    def test_hebrew_text_survives(self):
        converted = bibItem(["יוסף כוהן."]).convert()
        assert "יוסף" in converted

    @pytest.mark.xfail(
        strict=True,
        reason="B9: the math-mode branch tests `word[indx] == 'DOLLAR SIGN'`, "
               "comparing a single character to an 11-character string, so it can "
               "never fire. '$' falls through to the literal branch and the "
               "expression is mangled. See the consolidated findings issue.",
    )
    def test_math_mode_is_preserved(self):
        r"""The class docstring's own example is $4\frac{m}{s}$.

        It currently converts to $$4\$\L{frac}{\L{m}}{\L{s}}$ -- not valid
        LaTeX. Whatever the Hebrew wrapping does elsewhere, an existing math
        expression has to come out intact.
        """
        converted = bibItem([r"האם $4\frac{m}{s}$ נכון"]).convert()
        assert r"$4\frac{m}{s}$" in converted


@pytest.mark.unit
class TestBibtexFile:
    BIB = (
        "\\begin{thebibliography}{1}\n"
        "\\bibitem{one}\n"
        "John Smith.\n"
        "\\bibitem{two}\n"
        "Jane Doe.\n"
        "\\end{thebibliography}\n"
    )

    def test_reads_from_a_file_object(self):
        parsed = bibtexFile(io.StringIO(self.BIB))
        assert parsed.items  # non-empty

    def test_reads_from_a_path(self, tmp_path):
        path = tmp_path / "refs.bbl"
        path.write_text(self.BIB, encoding="utf-8")
        parsed = bibtexFile(str(path))
        assert parsed.items

    def test_missing_path_raises_valueerror(self, tmp_path):
        with pytest.raises(ValueError, match="does not exist"):
            bibtexFile(str(tmp_path / "absent.bbl"))

    def test_each_real_item_is_captured(self):
        parsed = bibtexFile(io.StringIO(self.BIB))
        joined = [ "".join(item._item) for item in parsed.items ]
        assert any("\\bibitem{one}" in text for text in joined)
        assert any("\\bibitem{two}" in text for text in joined)

    def test_convert_keeps_the_surrounding_environment(self):
        converted = bibtexFile(io.StringIO(self.BIB)).convert()
        assert converted.startswith("\\begin{thebibliography}")
        assert converted.rstrip().endswith("\\end{thebibliography}")

    def test_line_breaks_survive_conversion(self):
        """readlines() keeps the trailing newlines, so joining must not eat them."""
        converted = bibtexFile(io.StringIO(self.BIB)).convert()
        assert "John Smith." in converted
        assert converted.count("\n") > 2

    @pytest.mark.xfail(
        strict=True,
        reason="B10: _parseBibItems appends the accumulator on every \\bibitem "
               "line including the first, when it is still empty, so items[0] is "
               "always an empty bibItem. See the consolidated findings issue.",
    )
    def test_item_count_matches_the_bibitem_count(self):
        """Two \\bibitem lines must yield two items, not three."""
        parsed = bibtexFile(io.StringIO(self.BIB))
        assert len(parsed.items) == 2

    @pytest.mark.xfail(
        strict=True,
        reason="B11: an empty file indexes thefileStr[0] and raises IndexError "
               "instead of reporting an empty bibliography. "
               "See the consolidated findings issue.",
    )
    def test_empty_file_reports_clearly(self):
        with pytest.raises(ValueError):
            bibtexFile(io.StringIO(""))
