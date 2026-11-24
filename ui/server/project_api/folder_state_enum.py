"""What a project's files folder holds, as the manage list shows it."""
from enum import Enum


class FolderStateEnum(str, Enum):
    """The state of one project's files folder."""

    NO_FOLDER = "noFolder"  # No files folder is set.
    MISSING = "missing"  # A folder is set, but nothing is at that path.
    NOT_A_FOLDER = "notAFolder"  # Something is at that path, but it is a file.
    EMPTY = "empty"  # The folder holds nothing.
    ONLY_CASE_CONFIG = "onlyCaseConfig"  # The folder holds only caseConfiguration.json.
    HAS_FILES = "hasFiles"  # The folder holds other files or subfolders.
    UNREADABLE = "unreadable"  # The folder is there but cannot be listed.
