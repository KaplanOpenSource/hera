import { ListAlt } from "@mui/icons-material";
import { useState } from "react";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { ManageProjectsDialog } from "./ManageProjectsDialog";

// Icon button next to the project actions menu; opens the manage projects dialog.
export const ManageProjectsButton = () => {
  const [open, setOpen] = useState(false);

  const close = () => {
    setOpen(false);
  }

  return (<>
    <ButtonTooltip
      title='Manage projects'
      onClick={() => setOpen(true)}
    >
      <ListAlt />
    </ButtonTooltip>
    {open && (
      <ManageProjectsDialog
        onClose={close}
      />
    )}
  </>);
}
