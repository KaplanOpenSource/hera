"""_LazyModule: a proxy that imports its target on first attribute access.

The point of the class is that `import hera` must not drag in pint, pandas,
shapely and geopandas.  These tests check the two halves of that promise:
nothing is imported until something is asked for, and once it is, the proxy is
indistinguishable from the module.
"""
import sys

import pytest

from hera.utils.lazy import _LazyModule


@pytest.mark.unit
class TestDeferredImport:
    def test_construction_does_not_import(self):
        """The whole purpose: constructing the proxy must be free."""
        proxy = _LazyModule("base64")
        assert object.__getattribute__(proxy, "_mod") is None

    def test_repr_says_it_is_unloaded(self):
        proxy = _LazyModule("base64")
        assert "not yet loaded" in repr(proxy)
        assert "base64" in repr(proxy)

    def test_attribute_access_triggers_the_import(self):
        proxy = _LazyModule("base64")
        assert proxy.b64encode(b"hera") == b"aGVyYQ=="
        assert object.__getattribute__(proxy, "_mod") is sys.modules["base64"]

    def test_repr_becomes_the_module_repr_once_loaded(self):
        import base64

        proxy = _LazyModule("base64")
        proxy.b64encode  # force the load
        assert repr(proxy) == repr(base64)

    def test_the_module_is_loaded_only_once(self):
        proxy = _LazyModule("base64")
        proxy.b64encode
        first = object.__getattribute__(proxy, "_mod")
        proxy.b64decode
        assert object.__getattribute__(proxy, "_mod") is first


@pytest.mark.unit
class TestProxyTransparency:
    def test_attributes_match_the_real_module(self):
        import base64

        proxy = _LazyModule("base64")
        assert proxy.b64encode is base64.b64encode

    def test_missing_attribute_raises_attribute_error(self):
        proxy = _LazyModule("base64")
        with pytest.raises(AttributeError):
            proxy.no_such_function

    def test_a_missing_module_raises_on_first_access_not_construction(self):
        """Failure must be deferred too, otherwise the laziness is pointless."""
        proxy = _LazyModule("hera_module_that_does_not_exist")
        with pytest.raises(ModuleNotFoundError):
            proxy.anything

    def test_slots_prevent_attribute_assignment(self):
        """__slots__ is declared; setting a new attribute must fail.

        Without this the proxy would silently accept writes that never reach
        the underlying module.
        """
        proxy = _LazyModule("base64")
        with pytest.raises(AttributeError):
            object.__setattr__(proxy, "unexpected", 1)
