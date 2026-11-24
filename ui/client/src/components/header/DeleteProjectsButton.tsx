import { Button } from "@mui/material";
import { useState } from "react";
import { projectApi } from "../../api/api";
import { FolderStateEnum } from "../../shared/FolderStateEnum";
import { pushError, pushInfo } from "../../io/snackbar";
import { ProjectEntire, ProjectFolder, ProjectName } from "../../shared/types";
import { useProjectListStore } from "../../stores/useProjectListStore";
import { NO_PROJECT, useProjectStore } from "../../stores/useProjectStore";
import { DeleteProjectsConfirm } from "./DeleteProjectsConfirm";
import { ExportStep } from "./ExportStep";
import { exportProjectsToZip } from "./exportProjectsToZip";
import { projectsText } from "./projectsText";

// Saves the marked projects in projects.zip, then deletes them one by one.
export const DeleteProjectsButton = ({
  projects,
  onDeleted,
}: {
  projects: ProjectFolder[],
  onDeleted: (names: string[]) => void,
}) => {
  const [confirming, setConfirming] = useState(false);
  const [running, setRunning] = useState(false);
  const [step, setStep] = useState<ExportStep | null>(null);
  const { currProjectName, selectProject, setCurrentProject } = useProjectStore();
  const { setProjectNames } = useProjectListStore();

  // Shows the list left after the deletes, and leaves the current project if it is gone.
  const applyAnswer = (data: any, deleted: string[]) => {
    setProjectNames(data.projects as ProjectName[]);
    if (!deleted.includes(currProjectName)) {
      return;
    }
    const next = data.project as ProjectEntire | null;
    setCurrentProject(next);
    if (next) {
      selectProject(next.name);
    } else {
      selectProject(NO_PROJECT);
    }
  }

  const deleteAll = async () => {
    setRunning(true);
    const names: string[] = [];
    for (const project of projects) {
      names.push(project.name);
    }
    const failedExport = await exportProjectsToZip(names, setStep);
    if (failedExport) {
      pushError(`Delete stopped: ${failedExport} export failed.`);
      setRunning(false);
      setConfirming(false);
      return;
    }
    const deleted: string[] = [];
    let lastAnswer: any = null;
    for (const project of projects) {
      const name = project.name;
      // Only a folder known to be empty goes too; a folder with files is never deleted.
      const deleteFiles = project.folderState === FolderStateEnum.Empty
        || project.folderState === FolderStateEnum.OnlyCaseConfig;
      const answer = await projectApi.delete(name, deleteFiles);
      if (!answer) {
        pushError(`Delete stopped: ${name} delete failed. Deleted ${deleted.length} of ${projectsText(names.length)}.`);
        break;
      }
      deleted.push(name);
      lastAnswer = answer;
    }
    if (lastAnswer) {
      applyAnswer(lastAnswer, deleted);
    }
    if (deleted.length === names.length) {
      pushInfo(`Deleted ${projectsText(deleted.length)}.`);
    }
    onDeleted(deleted);
    setRunning(false);
    setConfirming(false);
  }

  return (<>
    <Button
      variant="outlined"
      color="error"
      disabled={projects.length === 0}
      onClick={() => setConfirming(true)}
      data-testid="manage-projects-delete"
    >
      {`Delete ${projectsText(projects.length)}...`}
    </Button>
    {confirming && (
      <DeleteProjectsConfirm
        projects={projects}
        running={running}
        step={step}
        onCancel={() => setConfirming(false)}
        onConfirm={deleteAll}
      />
    )}
  </>);
}
