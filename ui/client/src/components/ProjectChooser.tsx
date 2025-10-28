import { FormControl, InputLabel, Select, MenuItem, SelectChangeEvent, IconButton } from "@mui/material"
import { useEffect } from "react";
import { NO_PROJECT, useProjectStore } from "../stores/useProjectStore"
import { Add } from '@mui/icons-material';
import { execPython } from "../io/execPython";

export const ProjectChooser = ({ }) => {
  const { projectNames, currProjectName, selectProject } = useProjectStore();

  useEffect(() => {
    if (currProjectName === NO_PROJECT && projectNames.length > 0) {
      selectProject(projectNames[0].name);
    }
  }, [currProjectName, projectNames]);

  const addProject = () => {
    execPython(`
from types import SimpleNamespace
from hera.utils.data.CLI import project_create
project_create(SimpleNamespace(projectName="hello3", directory=None, loadRepositories=True, overwrite=False))
      `)
  }

  return (
    <>
      <FormControl size="small">
        <InputLabel id="demo-simple-select-label">Project</InputLabel>
        <Select
          labelId="demo-simple-select-label"
          id="demo-simple-select"
          value={currProjectName}
          label="Project"
          onChange={(event: SelectChangeEvent) => {
            selectProject(event.target.value as string);
          }}
        >
          {projectNames.map(({ name }) => (
            <MenuItem key={name} value={name}>
              {name}
            </MenuItem>
          ))}
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