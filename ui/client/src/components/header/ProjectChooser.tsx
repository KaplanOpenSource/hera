import { Autocomplete, TextField } from "@mui/material";
import { FolderOutlined } from "@mui/icons-material";
import { useEffect } from "react";
import { useNavigate } from "react-router-dom";
import { EMPTY_NAME_PROJECT, useProjectStore } from "../../stores/useProjectStore";
import { useProjectListStore } from "../../stores/useProjectListStore";
import { useServerUserStore } from "../../stores/useServerUserStore";

const displayName = (name: string) => name || EMPTY_NAME_PROJECT;
const storeName = (name: string) => name === EMPTY_NAME_PROJECT ? "" : name;

export const ProjectChooser = () => {
  const { projectNames } = useProjectListStore();
  const { currProjectName } = useProjectStore();
  const { username, loadServerUser } = useServerUserStore();
  const navigate = useNavigate();

  useEffect(() => {
    loadServerUser();
  }, [loadServerUser]);

  const options = projectNames.map(({ name }) => displayName(name));

  // The server user owns the projects, so name the list after them.
  let label = 'Projects';
  if (username) {
    label = `${username}'s projects`;
  }

  return (
    <Autocomplete
      size="small"
      value={displayName(currProjectName)}
      options={options}
      onChange={(_, value) => {
        if (value) {
          navigate('/' + encodeURIComponent(storeName(value)));
        }
      }}
      renderInput={(params) => (
        <TextField
          {...params}
          label={label}
          slotProps={{
            input: {
              ...params.InputProps,
              startAdornment: <FolderOutlined fontSize="small" sx={{ ml: 0.5, mr: 0.5, color: "text.secondary" }} />,
            },
          }}
        />
      )}
      disableClearable
      sx={{ minWidth: 260 }}
    />
  );
};
