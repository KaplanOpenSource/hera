import * as React from "react";
import {
  Dialog,
  DialogTitle,
  DialogContent,
  DialogContentText,
  DialogActions,
  Button,
  TextField,
} from "@mui/material";

export interface ConfirmOptions {
  title: string;
  message?: string;
  yesText?: string;
  noText?: string;
  /** If true, show a text input and return its value */
  requireText?: boolean;
  textLabel?: string;
  textPlaceholder?: string;
  textValidate?: (text: string) => boolean;
}

export interface ConfirmResult {
  confirmed: boolean;
  text?: string;
}

/**
 * Async confirm dialog hook (TypeScript)
 *
 * Usage:  
 *   const { confirm, ConfirmDialog } = useConfirm();  
 *   const result = await confirm({ title: 'Delete?', requireText: true });  
 */
export function useConfirm() {
  const [open, setOpen] = React.useState<boolean>(false);
  const [options, setOptions] = React.useState<ConfirmOptions>({ title: "" });
  const [text, setText] = React.useState<string>("");

  const resolverRef = React.useRef<((value: ConfirmResult) => void) | null>(null);

  const isValidated = options.textValidate ? options.textValidate(text) : true;

  const confirmOpen = React.useCallback((opts: ConfirmOptions): Promise<ConfirmResult> => {
    setOptions({
      yesText: "Yes",
      noText: "No",
      ...opts,
    });
    setText("");
    setOpen(true);

    return new Promise<ConfirmResult>((resolve) => {
      resolverRef.current = resolve;
    });
  }, []);

  const close = (confirmed: boolean) => {
    setOpen(false);
    resolverRef.current?.({ confirmed, text: options.requireText ? text : undefined });
    resolverRef.current = null;
  };

  const handleKeyDown = (event: React.KeyboardEvent) => {
    if (event.key === "Enter" && isValidated) {
      event.preventDefault();
      close(true);
    }
    if (event.key === "Escape") {
      event.preventDefault();
      close(false);
    }
  };

  const ConfirmDialog = (
    <Dialog
      open={open}
      onClose={() => close(false)}
      onKeyDown={handleKeyDown}
    >
      <DialogTitle>{options.title}</DialogTitle>
      <DialogContent>
        {options.message && (
          <DialogContentText sx={{ mb: options.requireText ? 2 : 0 }}>
            {options.message}
          </DialogContentText>
        )}

        {options.requireText && (
          <TextField
            autoFocus
            fullWidth
            label={options.textLabel ?? "Reason"}
            placeholder={options.textPlaceholder}
            value={text}
            onChange={(e) => setText(e.target.value)}
            onClick={(e) => e.stopPropagation()}
            sx={{ marginTop: 1 }}
            size='small'
            error={!isValidated}
          />
        )}
      </DialogContent>
      <DialogActions>
        <Button
          size='small'
          onClick={(e) => {
            e.stopPropagation();
            close(false);
          }}
        >
          {options.noText}
        </Button>
        <Button
          size='small'
          variant="contained"
          onClick={(e) => {
            e.stopPropagation();
            close(true);
          }}
          disabled={!isValidated}
        >
          {options.yesText}
        </Button>
      </DialogActions>
    </Dialog>
  );

  return { confirmOpen, ConfirmDialog } as const;
}

// ---------------- Example ----------------
// const { confirmOpen, ConfirmDialog } = useConfirm();
//
// const handleAction = async () => {
//   const result = await confirmOpen({
//     title: 'Delete item',
//     message: 'Please confirm deletion',
//     requireText: true,
//     textLabel: 'Type DELETE to confirm',
//   });
//
//   if (result.confirmed) {
//     console.log(result.text);
//   }
// };
//
// return (
//   <>
//     <Button onClick={handleAction}>Delete</Button>
//     {ConfirmDialog}
//   </>
// );
