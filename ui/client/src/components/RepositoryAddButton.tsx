import { LibraryAdd } from "@mui/icons-material";
import {
  Checkbox,
  FormControlLabel,
  FormGroup,
  Stack,
  TextField
} from "@mui/material";
import { ButtonTooltip } from "../elements/ButtonTooltip";
import { useDialog } from "../elements/useDialog";
import { execPython } from "../io/execPython";
import { useProjectStore } from "../stores/useProjectStore";
import { ProjectEntire } from "../shared/types";

export const RepositoryAddButton = ({ }) => {
  const { currProject, setCurrentProject } = useProjectStore();

  type RepoParams = {
    repositoryJson: string;
    baseDir: string;
    overwrite: boolean;
  };

  const { openDialog, DialogComponent } = useDialog<RepoParams>();

  const addRepo = async (params: RepoParams) => {
    const { data } = await execPython(`
import logging
import os
from hera.datalayer import All
from hera.utils.data.toolkit import dataToolkit
logger = logging.getLogger("hera.bin.repository_load")
dtk = dataToolkit()
dtk.loadAllDatasourcesInRepositoryJSONToProject(projectName='${currProject?.name}',
                                                repositoryJSON=${params.repositoryJson},
                                                basedir=os.path.abspath('${params.baseDir}'),
                                                overwrite=${params.overwrite ? 'True' : 'False'})
docs = All.getDocumentsAsDict('${currProject?.name}', with_id=True)
project = {"name": '${currProject?.name}', "documents": docs['documents']}
result = {"project": project}
      `);
    if (data) {
      setCurrentProject(data.project as ProjectEntire);
    }
  }

  return (<>
    <ButtonTooltip
      title={'Load Datasources To Project (Repository)'}
      onClick={async () => {
        const result = await openDialog({
          title: 'Load Datasources Repository To Project',
          initialValues: { repositoryJson: '', baseDir: '', overwrite: true },
          render: ({ values, setValues }) => (
            <Stack direction={'column'} spacing={2}>
              <TextField
                label="Repository Json (as string)"
                fullWidth
                value={values.repositoryJson}
                onChange={(e) => setValues({ ...values, repositoryJson: e.target.value })}
                size="small"
              />
              <TextField
                label="Base Directory"
                fullWidth
                value={values.baseDir}
                onChange={(e) => setValues({ ...values, baseDir: e.target.value })}
                size="small"
              />
              <FormGroup>
                <FormControlLabel
                  label="Overwrite"
                  control={<Checkbox
                    checked={values.overwrite}
                    onChange={(e) => setValues({ ...values, overwrite: e.target.checked })}
                  />}
                />
              </FormGroup>
            </Stack>
          ),
        });

        if (result.confirmed && result.values) {
          await addRepo(result.values)
        }
      }}
    >
      <LibraryAdd />
      {DialogComponent}
    </ButtonTooltip>
  </>)
}