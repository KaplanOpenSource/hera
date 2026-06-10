import { Settings } from "@mui/icons-material";
import { Button, createTheme, DialogActions, DialogContent, DialogTitle, Stack } from "@mui/material";
import { BooleanProperty } from "../../elements/BooleanProperty";
import { ButtonDialog } from "../../elements/ButtonDialog";
import { NumberProperty } from "../../elements/NumberProperty";
import { ReloadIntervalSlider } from "../../elements/ReloadIntervalSlider";
import { useViewSettingsStore } from "../../stores/useViewSettingsStore";

// Own (light) theme so the dialog isn't tinted by the dark app-header theme it opens from.
const dialogTheme = createTheme();

export const ProjectViewSettingsButton = ({ }) => {
  const { viewSettings, setViewSettings } = useViewSettingsStore();
  return (
    <ButtonDialog icon={<Settings />} title="Settings" dialogTheme={dialogTheme}>
      {(close) => (
        <>
          <DialogTitle>Settings</DialogTitle>
          <DialogContent>
            <Stack direction="column" spacing={1} sx={{ mt: 1 }}>
              <ReloadIntervalSlider
                value={viewSettings.reloadIntervalSeconds}
                setValue={(v) => setViewSettings({ ...viewSettings, reloadIntervalSeconds: v })}
              />
              <NumberProperty
                label="Max depth"
                value={viewSettings.maxDepth}
                min={0}
                max={25}
                setValue={(v) => setViewSettings({ ...viewSettings, maxDepth: v })}
              />
              <NumberProperty
                label="Min group size"
                value={viewSettings.minGroupSize}
                min={1}
                setValue={(v) => setViewSettings({ ...viewSettings, minGroupSize: v })}
              />
              <NumberProperty
                label="Max branches"
                value={viewSettings.maxBranches}
                min={1}
                setValue={(v) => setViewSettings({ ...viewSettings, maxBranches: v })}
              />
              <BooleanProperty
                label="First branch by header fields (toolkit, notebooks, type, etc...)"
                value={viewSettings.firstBranchHeadFields}
                setValue={v => setViewSettings({ ...viewSettings, firstBranchHeadFields: v })}
              />
              <BooleanProperty
                label="Show document preview"
                value={viewSettings.showDocumentPreview}
                setValue={v => setViewSettings({ ...viewSettings, showDocumentPreview: v })}
              />
            </Stack>
          </DialogContent>
          <DialogActions>
            <Button onClick={close}>Done</Button>
          </DialogActions>
        </>
      )}
    </ButtonDialog>
  )
}
