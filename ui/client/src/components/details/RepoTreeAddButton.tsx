import { DriveFolderUpload } from "@mui/icons-material";
import { BooleanProperty } from "../../elements/BooleanProperty";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { useDialog } from "../../elements/useDialog";
import { execPython } from "../../io/execPython";
import { ProjectEntire } from "../../shared/types";
import { useProjectStore } from "../../stores/useProjectStore";

export const RepoTreeAddButton = ({
  tree,
}: {
  tree: any,
}) => {
  type AddRepoArgs = {
    overwrite: boolean;
  };

  const { openDialog, DialogComponent } = useDialog<AddRepoArgs>();
  const { currProject, setCurrentProject } = useProjectStore();

  const baseDir = currProject?.documents.find(x => x.type === currProject.name + '__config__')?.desc.filesDirectory || '';

  const addRepo = async (params: AddRepoArgs) => {
    const { data } = await execPython(`
import logging
import os
from hera.datalayer import All
from hera.utils.data.toolkit import dataToolkit
logger = logging.getLogger("hera.bin.repository_load")
dtk = dataToolkit()
dtk.loadAllDatasourcesInRepositoryJSONToProject(projectName='${currProject?.name}',
                                                repositoryJSON=${JSON.stringify(tree)},
                                                basedir=os.path.abspath('${baseDir}'),
                                                overwrite=${params.overwrite ? 'True' : 'False'})
docs = All.getDocumentsAsDict('${currProject?.name}', with_id=True)
project = {"name": '${currProject?.name}', "documents": docs['documents']}
result = {"project": project}
        `);
    if (data) {
      setCurrentProject(data.project as ProjectEntire);
    }
  }


  return (
    <ButtonTooltip
      color="primary"
      title={'Add repository of data sources to project'}
      disabled={tree === undefined}
      onClick={async () => {
        const result = await openDialog({
          title: 'Add Repository of Datasources To Project',
          initialValues: { overwrite: true },
          render:
            ({ values, setValues }) => (
              <BooleanProperty
                label="Overwrite"
                value={values.overwrite}
                setValue={(v) => setValues({ ...values, overwrite: v })}
              />
            )
        });

        if (result.confirmed && result.values) {
          await addRepo(result.values)
        }
      }}
    >
      <DriveFolderUpload />
      {DialogComponent}
    </ButtonTooltip>
  )
}