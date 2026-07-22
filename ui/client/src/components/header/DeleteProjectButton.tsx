import { Delete } from "@mui/icons-material";
import { ThemeProvider } from "@mui/material";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { useConfirm } from "../../elements/useConfirm";
import { fetchPython } from "../../io/fetchPython";
import { ProjectEntire, ProjectName } from "../../shared/types";
import { useProjectStore } from "../../stores/useProjectStore";
import { useAppTheme } from "../../theme";

export const DeleteProjectButton = ({ }) => {
  const { confirmOpen, ConfirmDialog } = useConfirm()
  const { currProjectName, selectProject, setProjectNames, setCurrentProject } = useProjectStore();
  // Follow the app theme so the dialog isn't tinted by the dark header it opens from.
  const dialogTheme = useAppTheme();

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
  }

  return (<>
    <ButtonTooltip
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
      <Delete />
      <ThemeProvider theme={dialogTheme}>{ConfirmDialog}</ThemeProvider>
    </ButtonTooltip>
  </>)
}