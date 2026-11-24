import { projectApi } from "../../api/api";
import { CASE_CONFIGURATION } from "../../shared/constants";

// True only when the files folder is missing or holds nothing but the case configuration.
export const isFilesFolderEmpty = async (projectName: string) => {
  const data = await projectApi.files(projectName, '');
  if (!data) {
    return false;
  }
  for (const file of data.files) {
    if (file.name !== CASE_CONFIGURATION) {
      return false;
    }
  }
  return true;
}
