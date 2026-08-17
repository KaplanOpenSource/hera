import { Autocomplete, TextField } from "@mui/material";
import { FolderOutlined } from "@mui/icons-material";
import { useNavigate } from "react-router-dom";
import { EMPTY_NAME_PROJECT, useProjectStore } from "../../stores/useProjectStore";

const displayName = (name: string) => name || EMPTY_NAME_PROJECT;
const storeName = (name: string) => name === EMPTY_NAME_PROJECT ? "" : name;

export const ProjectChooser = () => {
  const { projectNames, currProjectName } = useProjectStore();
  const navigate = useNavigate();

  const options = projectNames.map(({ name }) => displayName(name));

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
          placeholder="Project"
          slotProps={{
            input: {
              ...params.InputProps,
              startAdornment: <FolderOutlined fontSize="small" sx={{ ml: 0.5, mr: 0.5, color: "text.secondary" }} />,
            },
          }}
        />
      )}
      disableClearable
      sx={{ minWidth: 200 }}
    />
  );
};