import { Settings } from "@mui/icons-material";
import { Stack, TextField } from "@mui/material";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { useDialog } from "../../elements/useDialog";
import { ViewSettingsType } from "./ProjectTreeView";

export const ProjectViewSettingsButton = ({
  viewSettings,
  setViewSettings,
}: {
  viewSettings: ViewSettingsType,
  setViewSettings: (newVal: ViewSettingsType) => void,
}) => {
  const { DialogComponent, openDialog } = useDialog();
  return (
    <ButtonTooltip
      title={'Change tree settings'}
      onClick={async () => {
        const result = await openDialog({
          title: 'Change tree settings',
          initialValues: { ...viewSettings },
          maxWidth: 'sm',
          render:
            ({ values, setValues }) => (
              <Stack direction={'column'} spacing={1}>
                <TextField
                  label="Max Depth"
                  value={values.maxDepth}
                  type='number'
                  variant='outlined'
                  slotProps={{ htmlInput: { min: 0, max: 25 } }}
                  size='small'
                  onChange={(e) => setValues({ ...values, maxDepth: e.target.value })}
                  onClick={e => e.stopPropagation()}
                />
                <TextField
                  label="Min Group Size"
                  value={values.minGroupSize}
                  type='number'
                  variant='outlined'
                  slotProps={{ htmlInput: { min: 1 } }}
                  size='small'
                  onChange={(e) => setValues({ ...values, minGroupSize: e.target.value })}
                  onClick={e => e.stopPropagation()}
                />
              </Stack>
            )
        });

        if (result.confirmed && result.values) {
          setViewSettings(result.values as ViewSettingsType);
        }
      }}
    >
      <Settings />
      {DialogComponent}
    </ButtonTooltip>
  )
}