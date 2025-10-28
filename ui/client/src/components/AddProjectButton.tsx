import { Add } from "@mui/icons-material"
import { IconButton } from "@mui/material"
import { execPython } from "../io/execPython";
import { fetchProjectsNames, fetchProjectDetails } from "../io/FetchProjects";
import { useProjectStore } from "../stores/useProjectStore";

export const AddProjectButton = ({ }) => {
  const { selectProject } = useProjectStore();

  const addProject = async () => {
    const name = prompt('New project name?');
    if (!name) {
      return;
    }
    const { problem } = await execPython(`
from types import SimpleNamespace
from hera.utils.data.CLI import project_create
project_create(SimpleNamespace(projectName="${name}", directory=None, loadRepositories=True, overwrite=False))
      `)
    if (problem) {
      return;
    }
    await fetchProjectsNames();
    await fetchProjectDetails(name);
    selectProject(name);
  }

  return (<>
    <IconButton
      onClick={() => addProject()}
    >
      <Add />
    </IconButton>
  </>)
}