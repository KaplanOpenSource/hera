import { zipSync } from "fflate";
import { projectApi } from "../../api/api";
import { saveFile } from "../../utils/saveFile";
import { ExportStep } from "./ExportStep";

// Answers the project's zip, or null when the export failed or threw.
const exportOne = async (name: string) => {
  try {
    return await projectApi.export(name);
  } catch (error) {
    console.error(`could not export ${name}`, error);
    return null;
  }
}

// Exports the projects one by one into one saved projects.zip; answers the first failed name, or null.
export const exportProjectsToZip = async (
  names: string[],
  onProgress: (step: ExportStep | null) => void,
) => {
  const zips: { [fileName: string]: Uint8Array } = {};
  let done = 0;
  for (const name of names) {
    onProgress({ done, total: names.length, current: name });
    const bytes = await exportOne(name);
    if (!bytes) {
      onProgress(null);
      return name;
    }
    zips[`${name}.zip`] = new Uint8Array(bytes);
    done += 1;
  }
  // The inner zips are already packed, so store them as they are.
  const zipped = zipSync(zips, { level: 0 });
  saveFile(zipped.slice().buffer, 'projects.zip');
  onProgress(null);
  return null;
}
