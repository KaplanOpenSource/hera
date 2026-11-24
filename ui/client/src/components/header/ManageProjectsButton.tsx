import { ListAlt } from "@mui/icons-material";
import { useState } from "react";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { ManageProjectsDialog } from "./ManageProjectsDialog";

// The project menu item that opens the manage projects dialog.
export const ManageProjectsButton = ({
  onClosed,
}: {
  onClosed?: () => void,
}) => {
  const [open, setOpen] = useState(false);

  const close = () => {
    setOpen(false);
    onClosed?.();
  }

  return (<>
    <ButtonTooltip
      button
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
