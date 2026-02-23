import { Settings } from "@mui/icons-material";
import { Stack, TextField } from "@mui/material";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { useDialog } from "../../elements/useDialog";
import { ViewSettingsType } from "./ProjectTreeView";
import { NumberProperty } from "../../elements/NumberProperty";

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
                <NumberProperty
                  label="Max Depth"
                  value={values.maxDepth}
                  min={0}
                  max={25}
                  setValue={v => setValues({ ...values, maxDepth: v })}
                />
                <NumberProperty
                  label="Min Group Size"
                  value={values.minGroupSize}
                  min={1}
                  setValue={v => setValues({ ...values, minGroupSize: v })}
                />
                <NumberProperty
                  label="Max Branches"
                  value={values.maxBranches}
                  min={1}
                  setValue={v => setValues({ ...values, maxBranches: v })}
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