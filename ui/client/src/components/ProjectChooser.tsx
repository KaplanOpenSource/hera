import { FormControl, InputLabel, Select, MenuItem, SelectChangeEvent, IconButton } from "@mui/material"
import { useEffect } from "react";
import { NO_PROJECT, useProjectStore } from "../stores/useProjectStore"
import { Add } from '@mui/icons-material';

export const ProjectChooser = ({ }) => {
  const { projectNames, currProjectName, selectProject } = useProjectStore();

  useEffect(() => {
    if (currProjectName === NO_PROJECT && projectNames.length > 0) {
      selectProject(projectNames[0].name);
    }
  }, [currProjectName, projectNames]);

  function addProject(): void {
    // throw new Error("Function not implemented.");
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