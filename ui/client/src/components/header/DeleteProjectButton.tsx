import { Delete } from "@mui/icons-material";
import { Checkbox, FormControlLabel, FormHelperText, Stack, TextField } from "@mui/material";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { useDialog } from "../../elements/useDialog";
import { fetchPython } from "../../io/fetchPython";
import { ProjectEntire, ProjectName } from "../../shared/types";
import { useProjectStore } from "../../stores/useProjectStore";

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
  const { currProjectName, selectProject, setProjectNames, setCurrentProject } = useProjectStore();

  const deleteProject = async (deleteFiles: boolean) => {
    const { data } = await fetchPython({
      results: ['projectNames', 'project'],
      label: `delete project ${currProjectName}`,
      code: `
import os, shutil
from hera.datalayer import All
from hera.datalayer.project import getProjectList, Project

deleteFiles = ${deleteFiles ? 'True' : 'False'}
filesDir = os.path.abspath(Project(projectName='${currProjectName}').filesDirectory or '')

docs = All.getDocumentsAsDict('${currProjectName}', with_id=True)['documents']
docs = sorted(docs, key=lambda d: d['type'] == '${currProjectName}__config__')
for d in docs:
    res = d.get('resource')
    # Only remove on-disk files that live inside the project's files directory.
    if deleteFiles and filesDir and isinstance(res, str) and res:
        resAbs = os.path.abspath(res)
        if (resAbs == filesDir or resAbs.startswith(filesDir + os.sep)) and os.path.exists(resAbs):
            shutil.rmtree(resAbs) if os.path.isdir(resAbs) else os.remove(resAbs)
    All.deleteDocumentByID(d['_id']['$oid'])
if deleteFiles and filesDir and os.path.isdir(filesDir):
    # Remove the caseConfiguration.json written at project creation, then drop the folder if it is now empty.
    caseConfig = os.path.join(filesDir, 'caseConfiguration.json')
    if os.path.exists(caseConfig):
        os.remove(caseConfig)
    if not os.listdir(filesDir):
        os.rmdir(filesDir)
projectNames = [{"name": proj} for proj in getProjectList() if proj != '${currProjectName}']

docs = All.getDocumentsAsDict(projectNames[0]['name'], with_id=True)
project = {"name": projectNames[0]['name'], "documents": docs['documents']}
`,
    });
    const details = (data?.project) as ProjectEntire;
    const names = (data?.projectNames || []) as ProjectName[];
    setProjectNames(names);
    setCurrentProject(details);
    selectProject(details.name);
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
