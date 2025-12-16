import {
  Checkbox,
  FormControlLabel,
  FormGroup,
  Stack,
  TextField,
  Typography,
} from "@mui/material"
import { useState } from "react";

export type RepoAddParams = {
  repositoryJson: string;
  baseDir: string;
  overwrite: boolean;
};

const parseJson = (str: string) => {
  try {
    return JSON.parse(str);
  } catch {
    return {}
  }
}
export const RepoAddEditor = ({
  values,
  setValues,
}: {
  values: RepoAddParams,
  setValues: (newValues: RepoAddParams) => void,
}) => {
  const [sourcesUsed, setSourcesUsed] = useState<string[]>([]);
  const jsonobj = parseJson(values.repositoryJson);
  const dataSources = Object.keys(jsonobj);
  return (
    <Stack direction={'column'} spacing={2}>
      <Stack direction={'row'} spacing={2}>
        <TextField
          label="Repository Json (as string)"
          value={values.repositoryJson}
          onClick={e => e.stopPropagation()}
          onChange={(e) => setValues({ ...values, repositoryJson: e.target.value })}
          size="small"
          rows={10}
          multiline={true}
        />
        <Stack>
          <Typography>
            DataSources
          </Typography>
          {dataSources.map(d => (
            <FormGroup key={d}>
              <FormControlLabel
                label={d}
                control={<Checkbox
                  checked={values.overwrite}
                  onChange={(e) => setValues({ ...values, overwrite: e.target.checked })}
                />}
              />
            </FormGroup>
          ))}
        </Stack>
      </Stack>
      <TextField
        label="Base Directory"
        fullWidth
        value={values.baseDir}
        onClick={e => e.stopPropagation()}
        onChange={(e) => setValues({ ...values, baseDir: e.target.value })}
        size="small"
      />
      <FormGroup>
        <FormControlLabel
          label="Overwrite"
          control={<Checkbox
            checked={values.overwrite}
            onChange={(e) => setValues({ ...values, overwrite: e.target.checked })}
          />}
        />
      </FormGroup>
    </Stack>
  )
}