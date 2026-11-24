import { FolderStateEnum } from "../../shared/FolderStateEnum";
import { ProjectFolder } from "../../shared/types";

// Nothing in it: no documents, and its folder is not set, missing, empty, or holds only the case config.
export const isEmptyProject = (project: ProjectFolder) => {
  const state = project.folderState;
  const noFiles = state === FolderStateEnum.NoFolder
    || state === FolderStateEnum.Missing
    || state === FolderStateEnum.Empty
    || state === FolderStateEnum.OnlyCaseConfig;
  return noFiles && project.documentCount === 0;
}
