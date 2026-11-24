import { Button } from "@mui/material";
import { zipSync } from "fflate";
import { useState } from "react";
import { projectApi } from "../../api/api";
import { pushError, pushInfo } from "../../io/snackbar";
import { saveFile } from "../../utils/saveFile";

// Answers the project's zip, or null when the export failed or threw.
const exportOne = async (name: string) => {
  try {
    return await projectApi.export(name);
  } catch (error) {
    console.error(`could not export ${name}`, error);
    return null;
  }
}

const projectsText = (count: number) => {
  if (count === 1) {
    return '1 project';
  }
  return `${count} projects`;
}

// Exports the marked projects one by one and saves them all in one projects.zip.
export const ExportProjectsButton = ({
  names,
}: {
  names: string[],
}) => {
  const [busy, setBusy] = useState(false);

  const exportAll = async () => {
    setBusy(true);
    const zips: { [fileName: string]: Uint8Array } = {};
    for (const name of names) {
      const bytes = await exportOne(name);
      if (!bytes) {
        pushError(`Export stopped: ${name} failed.`);
        setBusy(false);
        return;
      }
      zips[`${name}.zip`] = new Uint8Array(bytes);
    }
    // The inner zips are already packed, so store them as they are.
    const zipped = zipSync(zips, { level: 0 });
    saveFile(zipped.slice().buffer, 'projects.zip');
    pushInfo(`Exported ${projectsText(names.length)}.`);
    setBusy(false);
  }

  let label = `Export ${projectsText(names.length)}`;
  if (busy) {
    label = 'Exporting...';
  }
  const disabled = busy || names.length === 0;

  return (
    <Button
      variant="contained"
      disabled={disabled}
      onClick={exportAll}
      data-testid="manage-projects-export"
    >
      {label}
    </Button>
  );
}
