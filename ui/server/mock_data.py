from typing import List

from models import Project

MOCK_PROJECTS: List[Project] = [
    Project(
        id="p-1001",
        name="Apollo",
        documents={"sim": 12, "measure": 7, "cache": 3},
        toolkitCount=4,
    ),
    Project(
        id="p-1002",
        name="Hermes",
        documents={"sim": 5, "measure": 2, "cache": 9},
        toolkitCount=2,
    ),
    Project(
        id="p-1003",
        name="Athena",
        documents={"sim": 20, "measure": 15, "cache": 11},
        toolkitCount=6,
    ),
]
