import { LibraryAdd } from "@mui/icons-material"
import { ButtonTooltip } from "../elements/ButtonTooltip"
import { useDialog } from "../elements/useDialog";
import { Stack, TextField } from "@mui/material";

export const RepositoryAddButton = ({ }) => {
  const { openDialog, DialogComponent } = useDialog<{
    reason: string;
    count: number;
  }>();

  const handleClick = async () => {
    const result = await openDialog({
      title: 'Confirm action',
      initialValues: { reason: '', count: 1 },
      render: ({ values, setValue }) => (
        <Stack direction={'column'} spacing={1}>
          <TextField
            label="Reason"
            fullWidth
            value={values.reason}
            onChange={(e) => setValue('reason', e.target.value)}
          />
          <TextField
            type="number"
            label="Count"
            value={values.count}
            onChange={(e) => setValue('count', Number(e.target.value))}
          />
        </Stack>
      ),
    });

    if (result.confirmed) {
      console.log(result.values);
    }
  };

  return (<>
    <ButtonTooltip
      title={'Load Datasources To Project (Repository)'}
      onClick={handleClick}
    >
      <LibraryAdd />
      {DialogComponent}
    </ButtonTooltip>
  </>)
}