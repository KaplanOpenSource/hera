import { Keyboard } from "@mui/icons-material";
import {
  Box,
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  TextField,
} from "@mui/material";
import { useState } from "react";
import { ButtonTooltip } from "../../elements/ButtonTooltip";

export const EditAsJsonButton = ({
  data,
  setData,
}: {
  data: any;
  setData: (newVal: any) => void;
}) => {
  const [open, setOpen] = useState(false);
  const [jsonText, setJsonText] = useState("");
  const [isValid, setIsValid] = useState(true);
  const [isDirty, setIsDirty] = useState(false);
  const [confirmingExit, setConfirmingExit] = useState(false);

  const handleOpen = () => {
    const text = JSON.stringify(data, null, 2);
    setJsonText(text);
    setIsValid(true);
    setIsDirty(false);
    setConfirmingExit(false);
    setOpen(true);
  };

  const handleClose = () => {
    setOpen(false);
  };

  const handleChange = (e: React.ChangeEvent<HTMLTextAreaElement>) => {
    e.stopPropagation();
    const text = e.target.value;
    setJsonText(text);
    setConfirmingExit(false);

    let valid = true;
    try {
      JSON.parse(text);
    } catch {
      valid = false;
    }
    setIsValid(valid);
    setIsDirty(text !== JSON.stringify(data, null, 2));
  };

  const handleConfirm = () => {
    try {
      setData(JSON.parse(jsonText));
    } catch { }
    handleClose();
  };

  const handleCancel = () => {
    if (isDirty && !confirmingExit) {
      setConfirmingExit(true);
      return;
    }
    handleClose();
  };

  return (
    <>
      <ButtonTooltip title="Edit as Json Text" onClick={handleOpen}>
        <Keyboard />
      </ButtonTooltip>

      <Dialog
        open={open}
        onClose={confirmingExit ? undefined : handleCancel}
        fullWidth
        maxWidth="xl"
        onClick={(e) => e.stopPropagation()}
      >
        <DialogTitle>Edit as Json Text</DialogTitle>
        <DialogContent>
          <Box sx={{ marginTop: 1 }}>
            <TextField
              label="Json"
              fullWidth
              multiline
              minRows={4}
              maxRows={20}
              value={jsonText}
              onChange={handleChange}
              error={!isValid}
              helperText={isValid ? "" : "Invalid Json"}
              onClick={(e) => e.stopPropagation()}
            />
          </Box>
        </DialogContent>
        <DialogActions>
          <Button
            variant="outlined"
            onClick={handleCancel}
            color={confirmingExit ? "error" : "inherit"}
          >
            {confirmingExit ? "Are you sure?" : "No"}
          </Button>
          <Button
            variant="contained"
            disabled={!isValid || !isDirty}
            onClick={handleConfirm}
          >
            Yes
          </Button>
        </DialogActions>
      </Dialog>
    </>
  );
};