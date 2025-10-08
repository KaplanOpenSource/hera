import { FormControl, InputLabel, Select, MenuItem, SelectChangeEvent } from "@mui/material"
import { useEffect, useState } from "react";
import { useProjectStore } from "../stores/useProjectStore"

const NO_PROJECT = "* NONE *";
export const ProjectChooser = ({ }) => {
  const { projectNames } = useProjectStore();
  const [projectId, setProjectId] = useState<string>(NO_PROJECT);

  useEffect(() => {
    if (projectId === NO_PROJECT && projectNames.length > 0) {
      setProjectId(projectNames[0].id);
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
          setProjectId(event.target.value as string);
        }}
      >
        {projectNames.map(({ id, name }) => (
          <MenuItem key={id} value={id}>
            {name}
          </MenuItem>
        ))}
      </Select>
    </FormControl>
  )
}