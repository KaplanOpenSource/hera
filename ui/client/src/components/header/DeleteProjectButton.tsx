import { Delete } from "@mui/icons-material";
import { Stack } from "@mui/material";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { useConfirm } from "../../elements/useConfirm";
import { fetchPython } from "../../io/fetchPython";
import { ProjectEntire, ProjectName } from "../../shared/types";
import { useProjectStore } from "../../stores/useProjectStore";

export const DeleteProjectButton = ({
  onDeleted,
}: {
  onDeleted?: () => void,
}) => {
  const { confirmOpen, ConfirmDialog } = useConfirm()
  const { currProjectName, selectProject, setProjectNames, setCurrentProject } = useProjectStore();

  const deleteProject = async () => {
    const { data } = await fetchPython({
      results: ['projectNames', 'project'],
      label: `delete project ${currProjectName}`,
      code: `
from hera.datalayer import All
from hera.datalayer.project import getProjectList

docs = All.getDocumentsAsDict('${currProjectName}', with_id=True)['documents']
docs = sorted(docs, key=lambda d: d['type'] == '${currProjectName}__config__')
for d in docs:
    All.deleteDocumentByID(d['_id']['$oid'])
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
        const { confirmed, text } = await confirmOpen({
          title: `Are you sure you want to delete ${currProjectName}?`,
          requireText: true,
          textLabel: 'Type project name to confirm',
          textPlaceholder: currProjectName,
          textValidate: (text) => text === currProjectName,
        });
        if (confirmed && text === currProjectName) {
          await deleteProject();
        }
      }}
    >
      <Stack direction={'row'} alignItems={'center'} spacing={0.5}>
        <Delete fontSize="small" />
      </Stack>
    </ButtonTooltip>
    {ConfirmDialog}
  </>)
}