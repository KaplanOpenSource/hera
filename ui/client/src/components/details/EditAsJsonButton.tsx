import { Check, Close, Keyboard } from "@mui/icons-material";
import {
  Box,
  Dialog,
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
    } catch {}
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
        <DialogTitle
          sx={{
            display: "flex",
            alignItems: "center",
            justifyContent: "space-between",
          }}
        >
          Edit as Json Text
          <Box sx={{ display: "flex", gap: 0.5 }}>
            <ButtonTooltip
              title={confirmingExit ? "Are you sure?" : "Cancel"}
              onClick={handleCancel}
              color={confirmingExit ? "error" : "default"}
            >
              <Close />
            </ButtonTooltip>
            <ButtonTooltip
              title="Confirm"
              onClick={handleConfirm}
              color="primary"
              disabled={!isValid || !isDirty}
            >
              <Check />
            </ButtonTooltip>
          </Box>
        </DialogTitle>
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
      </Dialog>
    </>
  );
};