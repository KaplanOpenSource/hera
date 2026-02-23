import { Autocomplete, TextField, Stack } from "@mui/material";
import { EMPTY_NAME_PROJECT, useProjectStore } from "../../stores/useProjectStore";
import { AddProjectButton } from "./AddProjectButton";
import { DeleteProjectButton } from "./DeleteProjectButton";

const displayName = (name: string) => name || EMPTY_NAME_PROJECT;
const storeName = (name: string) => name === EMPTY_NAME_PROJECT ? "" : name;

export const ProjectChooser = () => {
  const { projectNames, currProjectName, selectProject } = useProjectStore();

  const options = projectNames.map(({ name }) => displayName(name));

  return (
    <Stack direction="row" spacing={1} alignItems="center">
      <Autocomplete
        size="small"
        value={displayName(currProjectName)}
        options={options}
        onChange={(_, value) => {
          if (value) selectProject(storeName(value));
        }}
        renderInput={(params) => <TextField {...params} label="Project" />}
        disableClearable
        sx={{ minWidth: 200 }}
      />
      <AddProjectButton />
      <DeleteProjectButton />
    </Stack>
  );
};