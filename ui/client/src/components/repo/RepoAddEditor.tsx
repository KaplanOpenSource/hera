import {
  Checkbox,
  FormControlLabel,
  FormGroup,
  Stack,
  TextField,
} from "@mui/material"

export type RepoAddParams = {
  repositoryJson: string;
  baseDir: string;
  overwrite: boolean;
};

export const RepoAddEditor = ({
  values,
  setValues,
}: {
  values: RepoAddParams,
  setValues: (newValues: RepoAddParams) => void,
}) => {
  return (
    <Stack direction={'column'} spacing={2}>
      <TextField
        label="Repository Json (as string)"
        fullWidth
        value={values.repositoryJson}
        onClick={e => e.stopPropagation()}
        onChange={(e) => setValues({ ...values, repositoryJson: e.target.value })}
        size="small"
        rows={10}
        multiline={true}
      />
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