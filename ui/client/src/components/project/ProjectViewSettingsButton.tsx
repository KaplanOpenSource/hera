import { Settings } from "@mui/icons-material";
import { Button, DialogActions, DialogContent, DialogTitle, Stack, Typography } from "@mui/material";
import { BooleanProperty } from "../../elements/BooleanProperty";
import { ButtonDialog } from "../../elements/ButtonDialog";
import { NumberProperty } from "../../elements/NumberProperty";
import { ReloadIntervalSlider } from "../../elements/ReloadIntervalSlider";
import { ThemeModeSwitch } from "../../elements/ThemeModeSwitch";
import { ThemeMode, useViewSettingsStore } from "../../stores/useViewSettingsStore";
import { useAppTheme } from "../../theme";

export const ProjectViewSettingsButton = ({ }) => {
  const { viewSettings, setViewSettings } = useViewSettingsStore();
  // Follow the app theme so the dialog isn't tinted by the dark header it opens from.
  const dialogTheme = useAppTheme();
  return (
    <ButtonDialog icon={<Settings />} title="Settings" dialogTheme={dialogTheme}>
      {(close) => (
        <>
          <DialogTitle>Settings</DialogTitle>
          <DialogContent>
            <Stack direction="column" spacing={1} sx={{ mt: 1 }}>
              <Stack direction="row" spacing={1} alignItems="center">
                <Typography>Dark mode</Typography>
                <ThemeModeSwitch
                  mode={viewSettings.themeMode}
                  setMode={(mode: ThemeMode) => setViewSettings({ ...viewSettings, themeMode: mode })}
                />
              </Stack>
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
              <BooleanProperty
                label="Auto save workflow before running"
                value={viewSettings.alwaysSaveBeforeRun}
                setValue={v => setViewSettings({ ...viewSettings, alwaysSaveBeforeRun: v })}
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
