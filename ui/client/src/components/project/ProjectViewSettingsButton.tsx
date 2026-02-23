import { Settings } from "@mui/icons-material";
import { Button, DialogActions, DialogContent, DialogTitle, Stack } from "@mui/material";
import { ButtonDialog } from "../../elements/ButtonDialog";
import { NumberProperty } from "../../elements/NumberProperty";
import { ViewSettingsType } from "./ProjectTreeView";
import { BooleanProperty } from "../../elements/BooleanProperty";

export const ProjectViewSettingsButton = ({
  viewSettings,
  setViewSettings,
}: {
  viewSettings: ViewSettingsType;
  setViewSettings: (newVal: ViewSettingsType) => void;
}) => (
  <ButtonDialog icon={<Settings />} title="Change tree settings">
    {(close) => (
      <>
        <DialogTitle>Change tree settings</DialogTitle>
        <DialogContent>
          <Stack direction="column" spacing={1} sx={{ mt: 1 }}>
            <NumberProperty
              label="Max Depth"
              value={viewSettings.maxDepth}
              min={0}
              max={25}
              setValue={(v) => setViewSettings({ ...viewSettings, maxDepth: v })}
            />
            <NumberProperty
              label="Min Group Size"
              value={viewSettings.minGroupSize}
              min={1}
              setValue={(v) => setViewSettings({ ...viewSettings, minGroupSize: v })}
            />
            <NumberProperty
              label="Max Branches"
              value={viewSettings.maxBranches}
              min={1}
              setValue={(v) => setViewSettings({ ...viewSettings, maxBranches: v })}
            />
            <BooleanProperty
              label="Show Empty Toolkits"
              value={viewSettings.showEmptyToolkits}
              setValue={v => setViewSettings({ ...viewSettings, showEmptyToolkits: v })}
            />
            <BooleanProperty
              label="Show Document Preview"
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
