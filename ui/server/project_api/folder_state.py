"""Find the state of one files folder."""
import os

from project_api.folder_state_enum import FolderStateEnum
from project_api.project_files import ProjectFiles


def folder_state(folder: str) -> FolderStateEnum:
    """Look at the folder's top level only, and say what it holds."""
    if not folder:
        return FolderStateEnum.NO_FOLDER
    if not os.path.exists(folder):
        return FolderStateEnum.MISSING
    if not os.path.isdir(folder):
        return FolderStateEnum.NOT_A_FOLDER
    try:
        names = os.listdir(folder)
    # A folder we may not list is shown as such, not as an error.
    except OSError:
        return FolderStateEnum.UNREADABLE
    if not names:
        return FolderStateEnum.EMPTY
    if names == [ProjectFiles.CASE_CONFIG]:
        return FolderStateEnum.ONLY_CASE_CONFIG
    return FolderStateEnum.HAS_FILES
