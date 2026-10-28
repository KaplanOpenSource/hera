"""The toolkit documents, which tell the UI which toolkits exist."""
from __future__ import annotations

from typing import Any, Dict, List


def toolkit_documents() -> List[Dict[str, Any]]:
    """Every registered toolkit's document, as the client parses them."""
    from hera import toolkitHome
    return toolkitHome.getToolkitDocuments()
