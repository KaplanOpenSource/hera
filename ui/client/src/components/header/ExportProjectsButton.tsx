import { Button } from "@mui/material";
import { useState } from "react";
import { pushError, pushInfo } from "../../io/snackbar";
import { ExportProgressDialog } from "./ExportProgressDialog";
import { ExportStep } from "./ExportStep";
import { exportProjectsToZip } from "./exportProjectsToZip";
import { projectsText } from "./projectsText";

// Exports the marked projects one by one and saves them all in one projects.zip.
export const ExportProjectsButton = ({
  names,
}: {
  names: string[],
}) => {
  const [busy, setBusy] = useState(false);
  const [step, setStep] = useState<ExportStep | null>(null);

  const exportAll = async () => {
    setBusy(true);
    const failed = await exportProjectsToZip(names, setStep);
    if (failed) {
      pushError(`Export stopped: ${failed} failed.`);
    } else {
      pushInfo(`Exported ${projectsText(names.length)}.`);
    }
    setBusy(false);
  }

  return (<>
    <Button
      variant="contained"
      disabled={names.length === 0}
      onClick={exportAll}
      data-testid="manage-projects-export"
    >
      {`Export ${projectsText(names.length)}`}
    </Button>
    {busy && <ExportProgressDialog step={step} />}
  </>);
}
