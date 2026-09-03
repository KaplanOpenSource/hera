"""simulations/utils/inputForModelsCreation.py: the jinja2 template
renderer used to produce model input files.

B142: ``render`` guards against an unset template or params map by
printing a message -- and then falls through to
``return renderedTemplate``, a name only ever bound inside the ``else``
branch. The guard therefore raises ``UnboundLocalError`` instead of
returning anything, so the "not set yet" path is strictly worse than no
guard at all: the user gets the message *and* a crash.
"""
import pytest

from hera.simulations.utils.inputForModelsCreation import InputForModelsCreator


@pytest.fixture()
def templates_dir(tmp_path):
    (tmp_path / "case.template").write_text("velocity {{ speed }};\nheight {{ z }};\n")
    return tmp_path


@pytest.fixture()
def creator(templates_dir):
    return InputForModelsCreator(str(templates_dir))


@pytest.mark.unit
class TestSetup:
    def test_nothing_is_set_on_a_fresh_instance(self, creator):
        assert creator.templateName is None
        assert creator.paramsMap is None
        assert creator.renderedTemplate is None

    def test_the_template_name_round_trips(self, creator):
        creator.setTemplate("case.template")
        assert creator.templateName == "case.template"

    def test_the_params_map_round_trips(self, creator):
        creator.setParamsMap({"speed": 5})
        assert creator.paramsMap == {"speed": 5}


@pytest.mark.unit
class TestRender:
    def test_it_substitutes_the_params_into_the_template(self, creator):
        creator.setTemplate("case.template")
        creator.setParamsMap({"speed": 5, "z": 10})
        assert creator.render() == "velocity 5;\nheight 10;"  # jinja2 drops the trailing newline

    def test_the_result_is_kept_on_the_instance(self, creator):
        creator.setTemplate("case.template")
        creator.setParamsMap({"speed": 5, "z": 10})
        rendered = creator.render()
        assert creator.renderedTemplate == rendered

    def test_a_missing_placeholder_renders_empty_rather_than_raising(self, creator):
        """jinja2's default undefined is silent, so an incomplete params
        map yields a blank rather than an error -- worth pinning, because
        it means a typo'd key produces a quietly malformed input file."""
        creator.setTemplate("case.template")
        creator.setParamsMap({"speed": 5})
        assert creator.render() == "velocity 5;\nheight ;"

    def test_a_save_path_writes_the_rendered_text_to_disk(self, creator, tmp_path):
        creator.setTemplate("case.template")
        creator.setParamsMap({"speed": 5, "z": 10})
        target = tmp_path / "out" / "case"
        target.parent.mkdir()
        rendered = creator.render(savePath=str(target))
        assert target.read_text() == rendered

    def test_an_unknown_template_name_raises(self, creator):
        creator.setTemplate("nosuch.template")
        creator.setParamsMap({})
        from jinja2 import TemplateNotFound

        with pytest.raises(TemplateNotFound):
            creator.render()


@pytest.mark.unit
class TestRenderGuardIsBroken:
    """B142: see the module docstring."""

    @pytest.mark.xfail(
        strict=True,
        reason="B142: render's 'not set yet' guard prints a message and "
               "then returns a name bound only inside the else branch, so "
               "it raises UnboundLocalError instead of returning. See the "
               "consolidated findings issue.",
    )
    def test_rendering_with_nothing_set_should_not_crash(self, creator):
        creator.render()

    def test_rendering_with_nothing_set_currently_raises(self, creator):
        """Characterisation of B142."""
        with pytest.raises(UnboundLocalError, match="renderedTemplate"):
            creator.render()

    def test_it_still_prints_the_warning_before_crashing(self, creator, capsys):
        """Characterisation of B142: the user gets the message *and* the
        traceback."""
        with pytest.raises(UnboundLocalError):
            creator.render()
        assert "not set yet" in capsys.readouterr().out

    def test_a_template_without_a_params_map_hits_the_same_guard(self, creator):
        """Characterisation of B142: either half being unset is enough."""
        creator.setTemplate("case.template")
        with pytest.raises(UnboundLocalError):
            creator.render()

    def test_a_params_map_without_a_template_hits_the_same_guard(self, creator):
        """Characterisation of B142."""
        creator.setParamsMap({"speed": 5})
        with pytest.raises(UnboundLocalError):
            creator.render()
