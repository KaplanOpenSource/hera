import { FormControl, InputLabel, MenuItem, Select, SelectChangeEvent } from "@mui/material";
import { EMPTY_NAME_PROJECT, useProjectStore } from "../stores/useProjectStore";
import { AddProjectButton } from "./AddProjectButton";

export const ProjectChooser = ({ }) => {
  const { projectNames, currProjectName, selectProject } = useProjectStore();

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
      <AddProjectButton
      />
    </>
  )
}