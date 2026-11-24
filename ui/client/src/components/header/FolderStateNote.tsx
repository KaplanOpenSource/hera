import { Typography } from "@mui/material";
import { FolderStateEnum } from "../../shared/FolderStateEnum";

// The short note for a project's folder state; states that need no note show nothing.
export const FolderStateNote = ({
  state,
}: {
  state: FolderStateEnum,
}) => {
  return (<>
    {state === FolderStateEnum.Missing && (
      <Typography variant="caption" color="warning.main" data-testid="manage-projects-missing">
        folder missing
      </Typography>
    )}
    {state === FolderStateEnum.NotAFolder && (
      <Typography variant="caption" color="warning.main" data-testid="manage-projects-not-a-folder">
        not a folder
      </Typography>
    )}
    {state === FolderStateEnum.Empty && (
      <Typography variant="caption" color="warning.main" data-testid="manage-projects-empty">
        empty folder
      </Typography>
    )}
    {state === FolderStateEnum.OnlyCaseConfig && (
      <Typography variant="caption" color="warning.main" data-testid="manage-projects-empty">
        empty folder (case config only)
      </Typography>
    )}
    {state === FolderStateEnum.Unreadable && (
      <Typography variant="caption" color="warning.main" data-testid="manage-projects-unreadable">
        cannot read folder
      </Typography>
    )}
  </>);
}
