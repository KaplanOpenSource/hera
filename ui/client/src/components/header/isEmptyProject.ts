import { ProjectFolder } from "../../shared/types";
import { FolderCheckEnum } from "./FolderCheckEnum";

// Nothing in it: no documents, and its folder is empty, missing, or not set.
export const isEmptyProject = (project: ProjectFolder, folderCheck: FolderCheckEnum) => {
  const noFiles = folderCheck === FolderCheckEnum.Empty
    || folderCheck === FolderCheckEnum.Missing
    || folderCheck === FolderCheckEnum.NoFolder;
  return noFiles && project.documentCount === 0;
}
