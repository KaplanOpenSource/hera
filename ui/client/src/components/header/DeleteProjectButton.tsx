import { Delete } from "@mui/icons-material";
import { Checkbox, FormControlLabel, FormHelperText, Stack, TextField } from "@mui/material";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { useDialog } from "../../elements/useDialog";
import { postJson } from "../../io/entrypoint";
import { ProjectEntire, ProjectName } from "../../shared/types";
import { NO_PROJECT, useProjectStore } from "../../stores/useProjectStore";
import { useProjectListStore } from "../../stores/useProjectListStore";

interface DeleteProjectValues {
  confirmName: string;
  deleteFiles: boolean;
}

export const DeleteProjectButton = ({
  onDeleted,
}: {
  onDeleted?: () => void,
}) => {
  const { openDialog, DialogComponent } = useDialog<DeleteProjectValues>();
  const { currProjectName, selectProject, setCurrentProject } = useProjectStore();
  const { setProjectNames } = useProjectListStore();

  const deleteProject = async (deleteFiles: boolean) => {
    const data = await postJson('/project/delete', {
      projectName: currProjectName,
      deleteFiles,
    }, { label: `delete project ${currProjectName}` });
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
        const { confirmed, values } = await openDialog({
          title: `Are you sure you want to delete ${currProjectName}?`,
          maxWidth: 'sm',
          initialValues: { confirmName: '', deleteFiles: false },
          render: ({ values, setValues }) => (
            <Stack spacing={2}>
              <TextField
                autoFocus
                fullWidth
                size='small'
                label='Type project name to confirm'
                placeholder={currProjectName}
                value={values.confirmName}
                onChange={(e) => setValues({ ...values, confirmName: e.target.value })}
                onClick={(e) => e.stopPropagation()}
                error={values.confirmName !== '' && values.confirmName !== currProjectName}
              />
              <div>
                <FormControlLabel
                  control={
                    <Checkbox
                      checked={values.deleteFiles}
                      onChange={(e) => setValues({ ...values, deleteFiles: e.target.checked })}
                    />
                  }
                  label="Also delete the project's files from disk"
                />
                <FormHelperText>
                  Removes files under the project's files directory. External data referenced from elsewhere is left untouched.
                </FormHelperText>
              </div>
            </Stack>
          ),
        });
        if (confirmed && values?.confirmName === currProjectName) {
          await deleteProject(values.deleteFiles);
        }
      }}
    >
      <Stack direction={'row'} alignItems={'center'} spacing={0.5}>
        <Delete fontSize="small" />
      </Stack>
    </ButtonTooltip>
    {DialogComponent}
  </>)
}
