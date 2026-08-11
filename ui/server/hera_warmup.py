import logging
import threading

# Reuse uvicorn's logger so warmup lines match the rest of the server output.
logger = logging.getLogger("uvicorn.error")


class HeraWarmup:
    """Warms hera's slow imports off-thread so the server can serve index.html immediately (#1011)."""

    def __init__(self) -> None:
        self._ready = False

    def start(self) -> None:
        threading.Thread(target=self._warm, name="warm-hera", daemon=True).start()

    def _warm(self) -> None:
        logger.info("hera warmup starting")
        from hera import toolkitHome  # noqa: F401
        self._ready = True
        logger.info("hera warmup complete — server ready")

    @property
    def ready(self) -> bool:
        return self._ready
