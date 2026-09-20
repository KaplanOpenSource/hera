import time


def now_readable() -> str:
    """Local time as ``YYYY-MM-DD HH:MM:SS.mmm`` (millis so quick polls differ)."""
    now = time.time()
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(now)) + f".{int((now % 1) * 1000):03d}"
