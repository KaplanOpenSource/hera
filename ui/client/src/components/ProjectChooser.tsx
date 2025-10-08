import { FormControl, InputLabel, Select, MenuItem, SelectChangeEvent } from "@mui/material"
import { useEffect } from "react";
import { NO_PROJECT, useProjectStore } from "../stores/useProjectStore"

export const ProjectChooser = ({ }) => {
  const { projectNames, currProjectName: projectId, selectProject } = useProjectStore();

  useEffect(() => {
    if (projectId === NO_PROJECT && projectNames.length > 0) {
      selectProject(projectNames[0].name);
    }
  }, [projectId, projectNames]);

  return (
    <FormControl size="small">
      <InputLabel id="demo-simple-select-label">Project</InputLabel>
      <Select
        labelId="demo-simple-select-label"
        id="demo-simple-select"
        value={projectId}
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
  )
}