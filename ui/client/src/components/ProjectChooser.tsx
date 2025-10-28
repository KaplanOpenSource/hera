import { FormControl, InputLabel, Select, MenuItem, SelectChangeEvent, IconButton } from "@mui/material"
import { EMPTY_NAME_PROJECT, useProjectStore } from "../stores/useProjectStore"
import { Add } from '@mui/icons-material';
import { execPython } from "../io/execPython";
import { fetchProjectDetails, fetchProjectsNames } from "../io/FetchProjects";

export const ProjectChooser = ({ }) => {
  const { projectNames, currProjectName, selectProject } = useProjectStore();

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
  }

  const unempty = (s: string) => s === '' ? EMPTY_NAME_PROJECT : s;
  const reempty = (s: string) => s === EMPTY_NAME_PROJECT ? '' : s;

  return (
    <>
      <FormControl size="small">
        <InputLabel id="demo-simple-select-label">Project</InputLabel>
        <Select
          labelId="demo-simple-select-label"
          id="demo-simple-select"
          value={unempty(currProjectName)}
          label="Project"
          onChange={(event: SelectChangeEvent) => {
            selectProject(reempty(event.target.value as string));
          }}
        >
          {projectNames.map(({ name }) => {
            return (
              <MenuItem key={unempty(name)} value={unempty(name)}>
                {unempty(name)}
              </MenuItem>
            )
          })}
        </Select>
      </FormControl>
      <IconButton
        onClick={() => addProject()}
      >
        <Add />
      </IconButton>
    </>
  )
}