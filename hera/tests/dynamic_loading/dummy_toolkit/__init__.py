from typing import Any, Dict


class DummyToolkit:
    """
    Very small dynamic toolkit used only by tests.

    It can be instantiated with arbitrary kwargs and has a simple ping()
    method so we can verify it was loaded correctly.
    """
    def __init__(self, projectName: str = None, **kwargs: Any) -> None:
        self.projectName = projectName
        self.init_kwargs: Dict[str, Any] = dict(kwargs)

    def ping(self) -> str:
        return f"dummy-toolkit-ok:{self.projectName}"
