"""Helpers shared by the project tests: zips, the load route, notebook records."""
import asyncio
import io
import json
import zipfile

from routers import projects as routes


def _notebook(path):
    return {"_id": {"$oid": "nb"}, "type": "notebook", "resource": path}


def _zip_of(*entries):
    """A zip like hera's export, one JSON list of records per entry."""
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, mode="w") as zf:
        for i, records in enumerate(entries):
            zf.writestr(f"chunk_{i}", json.dumps(records))
    return buffer.getvalue()


class FakeRequest:
    """Only the raw body, which is all the load route reads."""

    def __init__(self, body):
        self._body = body

    async def body(self):
        return self._body


def _load(body):
    return asyncio.run(routes.load(FakeRequest(body)))
