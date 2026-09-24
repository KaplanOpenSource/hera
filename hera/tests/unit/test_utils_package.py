"""hera.utils lazy package surface, and the logging name helpers.

hera/utils/__init__.py exists so that `import hera` does not drag in pint,
pandas, shapely and geopandas.  Its docstring makes that promise; these tests
check both halves of it -- that every advertised name still resolves through
the PEP 562 __getattr__, and that a wrong name fails cleanly.
"""
import logging

import pytest

import hera.utils as heraUtils
from hera.utils.logging.helpers import (
    getClassLogger,
    get_classMethod_logger,
    get_logger,
    with_logger,
)


@pytest.mark.unit
class TestLazyPackageSurface:
    def test_every_advertised_name_resolves(self):
        """__all__ is the package's contract; a stale entry is a broken import.

        Checked as one assertion over the whole list so a single missing name
        names itself instead of hiding behind the first failure.
        """
        unresolved = [
            name for name in heraUtils.__all__ if not hasattr(heraUtils, name)
        ]
        assert unresolved == []

    def test_all_has_no_duplicates(self):
        assert len(heraUtils.__all__) == len(set(heraUtils.__all__))

    def test_star_import_names_are_importable(self):
        """`from hera.utils import *` must not raise on any advertised name."""
        namespace = {}
        exec("from hera.utils import *", namespace)  # noqa: S102 - deliberate
        for name in heraUtils.__all__:
            assert name in namespace, f"{name} is in __all__ but did not star-import"

    def test_unknown_name_raises_attribute_error(self):
        with pytest.raises(AttributeError, match="has no attribute"):
            heraUtils.no_such_helper

    def test_private_name_is_rejected_without_searching(self):
        """Underscore names short-circuit, so a typo cannot trigger imports."""
        with pytest.raises(AttributeError, match="has no attribute"):
            heraUtils._not_a_real_private

    def test_dir_includes_the_advertised_names(self):
        listed = set(heraUtils.__dir__())
        assert set(heraUtils.__all__) <= listed

    def test_dir_is_sorted_and_unique(self):
        listed = heraUtils.__dir__()
        assert listed == sorted(set(listed))

    def test_logging_helpers_are_eager(self):
        """The docstring keeps logging eager; it is needed during `import hera`."""
        assert heraUtils.get_logger is get_logger

    def test_a_resolved_name_is_cached_on_the_module(self):
        """__getattr__ writes into globals(), so the second lookup is direct."""
        first = heraUtils.toMeteorologicalAngle
        assert "toMeteorologicalAngle" in vars(heraUtils)
        assert heraUtils.toMeteorologicalAngle is first

    @pytest.mark.parametrize(
        "name, owning_submodule",
        [
            ("ureg", "unitHandler"),
            ("ConfigurationToJSON", "jsonutils"),
            ("dictToMongoQuery", "query"),
            ("toMeteorologicalAngle", "angle"),
            ("zip_items", "zipUtils"),
            ("standardize_polygon", "matplotlibCountour"),
        ],
    )
    def test_names_resolve_from_their_owning_submodule(self, name, owning_submodule):
        import importlib

        resolved = getattr(heraUtils, name)
        owner = importlib.import_module(f"hera.utils.{owning_submodule}")
        assert resolved is getattr(owner, name)


@pytest.mark.unit
class TestLoggerNaming:
    """Logger names must be module-qualified, so config can target them."""

    class Sample:
        pass

    def test_class_logger_is_module_qualified(self):
        logger = getClassLogger(self.Sample)
        assert logger.name.endswith("TestLoggerNaming.Sample")
        assert logger.name.startswith(self.Sample.__module__)

    def test_get_logger_from_an_instance(self):
        logger = get_logger(self.Sample())
        assert logger.name == getClassLogger(self.Sample).name

    def test_get_logger_with_an_explicit_name_ignores_the_instance(self):
        logger = get_logger(self.Sample(), "explicit.name")
        assert logger.name == "explicit.name"

    def test_class_method_logger_appends_the_method(self):
        logger = get_classMethod_logger(self.Sample(), "doThing")
        assert logger.name == f"{getClassLogger(self.Sample).name}.doThing"

    def test_class_method_logger_without_a_name_matches_the_class_logger(self):
        logger = get_classMethod_logger(self.Sample())
        assert logger.name == getClassLogger(self.Sample).name

    def test_the_same_name_returns_the_same_logger(self):
        """logging caches by name; two calls must not make two loggers."""
        assert get_logger(self.Sample()) is get_logger(self.Sample())

    def test_returns_real_logger_objects(self):
        assert isinstance(get_logger(self.Sample()), logging.Logger)


@pytest.mark.unit
class TestWithLogger:
    """Builds the (name, config, section) triple that initialize_logging consumes."""

    def test_returns_a_name_a_config_and_the_section(self):
        name, config, section = with_logger("hera.test")
        assert name == "hera.test"
        assert isinstance(config, dict)
        assert section == "loggers"

    def test_level_is_included_when_given(self):
        _, config, _ = with_logger("hera.test", level="DEBUG")
        assert config["level"] == "DEBUG"

    def test_handlers_are_included_when_given(self):
        _, config, _ = with_logger("hera.test", handlers=["console"])
        assert config["handlers"] == ["console"]

    def test_propagate_false_is_kept(self):
        """Only None means 'unspecified'; propagate=False is a real setting.

        A filter written as `if not value` instead of `if value is None` would
        silently drop this, which is why it gets its own test.
        """
        _, config, _ = with_logger("hera.test", propagate=False)
        assert config["propagate"] is False

    def test_unspecified_keys_are_omitted(self):
        """An absent key means 'inherit', which is not the same as a null value."""
        _, config, _ = with_logger("hera.test")
        assert config == {}

    @pytest.mark.xfail(
        strict=True,
        reason="B14: with_logger is annotated `-> (str, dict)` but returns a "
               "3-tuple (name, config, 'loggers'). Anyone unpacking it per the "
               "annotation gets ValueError: too many values to unpack. "
               "See the consolidated findings issue.",
    )
    def test_the_return_annotation_matches_the_arity(self):
        annotated = with_logger.__annotations__["return"]
        assert len(annotated) == len(with_logger("hera.test"))
