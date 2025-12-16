import { LibraryAdd } from "@mui/icons-material";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { useDialog } from "../../elements/useDialog";
import { execPython } from "../../io/execPython";
import { ProjectEntire } from "../../shared/types";
import { useProjectStore } from "../../stores/useProjectStore";
import { RepoAddEditor, RepoAddParams } from "./RepoAddEditor";

export const RepoAddButton = ({ }) => {
  const { currProject, setCurrentProject } = useProjectStore();

  const { openDialog, DialogComponent } = useDialog<RepoAddParams>();

  const addRepo = async (params: RepoAddParams) => {
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

  const baseDir = currProject?.documents.find(x => x.type === currProject.name + '__config__')?.desc.filesDirectory || '';

  return (<>
    <ButtonTooltip
      title={'Load Datasources To Project (Repository)'}
      onClick={async () => {
        const result = await openDialog({
          title: 'Load Datasources Repository To Project',
          initialValues: { repositoryJson: '', baseDir: baseDir, overwrite: true },
          render:
            ({ values, setValues }) => (
              <RepoAddEditor
                values={values}
                setValues={setValues}
              />
            )
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