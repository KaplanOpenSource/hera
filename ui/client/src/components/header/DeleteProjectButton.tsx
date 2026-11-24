import { Delete } from "@mui/icons-material";
import { Checkbox, FormControlLabel, FormHelperText, Radio, RadioGroup, Stack, TextField } from "@mui/material";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { useDialog } from "../../elements/useDialog";
import { projectApi } from "../../api/api";
import { ProjectEntire, ProjectName } from "../../shared/types";
import { NO_PROJECT, useProjectStore } from "../../stores/useProjectStore";
import { useProjectListStore } from "../../stores/useProjectListStore";
import { pushError } from "../../io/snackbar";
import { saveFile } from "../../utils/saveFile";
import { isFilesFolderEmpty } from "./isFilesFolderEmpty";

// How the user confirms the delete: by typing the name, or by saving a zip first.
enum ConfirmBy {
  Name = 'name',
  Export = 'export',
}

interface DeleteProjectValues {
  confirmName: string;
  deleteFiles: boolean;
  confirmBy: ConfirmBy;
}

export const DeleteProjectButton = ({
  onDeleted,
}: {
  onDeleted?: () => void,
}) => {
  const { openDialog, DialogComponent } = useDialog<DeleteProjectValues>();
  const { currProjectName, selectProject, setCurrentProject } = useProjectStore();
  const { setProjectNames } = useProjectListStore();

  // Saves the zip first; answers false when the export failed.
  const exportProject = async () => {
    const bytes = await projectApi.export(currProjectName);
    if (!bytes) {
      pushError('Export failed, so the project was not deleted. Try again.');
      return false;
    }
    saveFile(bytes, `${currProjectName}.zip`);
    return true;
  }

  const deleteProject = async (deleteFiles: boolean, exportFirst: boolean) => {
    if (exportFirst && !await exportProject()) {
      return;
    }
    const data = await projectApi.delete(currProjectName, deleteFiles);
    if (!data) {
      return;
    }
    setProjectNames(data.projects as ProjectName[]);
    const next = data.project as ProjectEntire | null;
    setCurrentProject(next);
    // Nothing is left to show when that was the last project.
    if (next) {
      selectProject(next.name);
    } else {
      selectProject(NO_PROJECT);
    }
    onDeleted?.();
  }

  return (<>
    <ButtonTooltip
      button
      color="error"
      title='Delete project'
      onClick={async () => {
        const filesEmpty = await isFilesFolderEmpty(currProjectName);
        const { confirmed, values } = await openDialog({
          title: `Are you sure you want to delete ${currProjectName}?`,
          maxWidth: 'sm',
          initialValues: { confirmName: '', deleteFiles: false, confirmBy: ConfirmBy.Name },
          render: ({ values, setValues }) => (
            <Stack spacing={2}>
              <RadioGroup
                value={values.confirmBy}
                onChange={(e) => setValues({ ...values, confirmBy: e.target.value as ConfirmBy })}
              >
                <FormControlLabel value={ConfirmBy.Name} control={<Radio />} label='Type the project name to confirm' />
                <TextField
                  autoFocus
                  fullWidth
                  size='small'
                  label='Project name'
                  placeholder={currProjectName}
                  disabled={values.confirmBy !== ConfirmBy.Name}
                  value={values.confirmName}
                  onChange={(e) => setValues({ ...values, confirmName: e.target.value })}
                  onClick={(e) => e.stopPropagation()}
                  error={values.confirmName !== '' && values.confirmName !== currProjectName}
                />
                <FormControlLabel
                  value={ConfirmBy.Export}
                  control={<Radio />}
                  disabled={values.deleteFiles}
                  label='Or export the project as a zip first'
                />
                <FormHelperText disabled={values.deleteFiles}>
                  You can restore the project from this zip later. Open the project menu, click 'Load project', and pick the zip.
                </FormHelperText>
                {values.deleteFiles && <FormHelperText>
                  Not possible when the files are deleted too.
                </FormHelperText>}
              </RadioGroup>
              <div>
                <FormControlLabel
                  control={
                    <Checkbox
                      checked={values.deleteFiles || filesEmpty}
                      disabled={filesEmpty}
                      onChange={(e) => setValues({ ...values, deleteFiles: e.target.checked, confirmBy: ConfirmBy.Name })}
                    />
                  }
                  label="Also delete the project's files from disk"
                />
                <FormHelperText>
                  Removes files under the project's files directory. External data referenced from elsewhere is left untouched.
                </FormHelperText>
                {filesEmpty && <FormHelperText>
                  The project's files folder is empty, so there is nothing to delete.
                </FormHelperText>}
              </div>
            </Stack>
          ),
        });
        if (!confirmed || !values) {
          return;
        }
        if (values.confirmBy === ConfirmBy.Name && values.confirmName !== currProjectName) {
          return;
        }
        // An empty folder is removed too, so no stale case configuration is left behind.
        await deleteProject(values.deleteFiles || filesEmpty, values.confirmBy === ConfirmBy.Export);
      }}
    >
      <Stack direction={'row'} alignItems={'center'} spacing={0.5}>
        <Delete fontSize="small" />
      </Stack>
    </ButtonTooltip>
    {DialogComponent}
  </>)
}
