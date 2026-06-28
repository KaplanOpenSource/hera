import importlib


class _LazyModule:
    """Proxy that imports a module the first time any attribute is accessed.

    Usage::

        numpy = _LazyModule("numpy")
        np    = numpy          # alias — same proxy object
        pd    = _LazyModule("pandas")

    The module is loaded once and cached; subsequent attribute accesses are
    direct attribute lookups on the real module.
    """
    __slots__ = ("_name", "_mod")

    def __init__(self, name):
        object.__setattr__(self, "_name", name)
        object.__setattr__(self, "_mod", None)

    def _load(self):
        mod = importlib.import_module(object.__getattribute__(self, "_name"))
        object.__setattr__(self, "_mod", mod)
        return mod

    def __getattr__(self, item):
        mod = object.__getattribute__(self, "_mod") or self._load()
        return getattr(mod, item)

    def __call__(self, *a, **kw):
        mod = object.__getattribute__(self, "_mod") or self._load()
        return mod(*a, **kw)

    def __repr__(self):
        mod = object.__getattribute__(self, "_mod")
        name = object.__getattribute__(self, "_name")
        if mod is None:
            return f"<_LazyModule '{name}' (not yet loaded)>"
        return repr(mod)
