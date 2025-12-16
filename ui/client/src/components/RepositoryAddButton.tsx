import { LibraryAdd } from "@mui/icons-material";
import {
  Checkbox,
  FormControlLabel,
  FormGroup,
  Stack,
  TextField
} from "@mui/material";
import { ButtonTooltip } from "../elements/ButtonTooltip";
import { useDialog } from "../elements/useDialog";

export const RepositoryAddButton = ({ }) => {
  const { openDialog, DialogComponent } = useDialog<{
    repositoryJson: string;
    baseDir: string;
    overwrite: boolean;
  }>();


  return (<>
    <ButtonTooltip
      title={'Load Datasources To Project (Repository)'}
      onClick={async () => {
        const result = await openDialog({
          title: 'Load Datasources Repository To Project',
          initialValues: { repositoryJson: '', baseDir: '', overwrite: true },
          render: ({ values, setValues }) => (
            <Stack direction={'column'} spacing={2}>
              <TextField
                label="Repository Json (as string)"
                fullWidth
                value={values.repositoryJson}
                onChange={(e) => setValues({ ...values, repositoryJson: e.target.value })}
                size="small"
              />
              <TextField
                label="Base Directory"
                fullWidth
                value={values.baseDir}
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
          ),
        });

        if (result.confirmed) {
          console.log(result.values);
        }
      }}
    >
      <LibraryAdd />
      {DialogComponent}
    </ButtonTooltip>
  </>)
}