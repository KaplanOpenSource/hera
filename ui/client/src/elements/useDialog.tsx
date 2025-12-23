import {
  Box,
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
} from "@mui/material";
import {
  ReactNode,
  useCallback,
  useRef,
  useState
} from "react";

export interface DialogOptions<TValues> {
  title: string;
  yesText?: string;
  noText?: string;
  /** Initial values for the dialog form */
  initialValues: TValues;
  /** Render function for dialog body */
  render: (args: {
    values: TValues;
    setValues: (newValues: TValues) => void;
  }) => ReactNode;
}

export interface DialogResult<TValues> {
  confirmed: boolean;
  values?: TValues;
}

/**
 * Generic async dialog hook
 *
 * Fully typed form values
 * Returns all values when confirmed
 */
export function useDialog<TValues extends Record<string, any>>() {
  const [open, setOpen] = useState(false);
  const [options, setOptions] = useState<DialogOptions<TValues> | null>(null);
  const [values, setValues] = useState<TValues | null>(null);

  const resolverRef = useRef<((value: DialogResult<TValues>) => void) | null>(null);

  const openDialog = useCallback(
    (opts: DialogOptions<TValues>): Promise<DialogResult<TValues>> => {
      setOptions({
        yesText: "Yes",
        noText: "No",
        ...opts,
      });
      setValues(opts.initialValues);
      setOpen(true);

      return new Promise((resolve) => {
        resolverRef.current = resolve;
      });
    },
    []
  );

  const close = (confirmed: boolean) => {
    setOpen(false);

    queueMicrotask(() => {
      resolverRef.current?.({
        confirmed,
        values: confirmed ? values ?? undefined : undefined,
      });
      resolverRef.current = null;
    });
  };

  const DialogComponent = options ? (
    <Dialog open={open} onClose={() => close(false)} fullWidth maxWidth={"xl"}>
      <DialogTitle>{options.title}</DialogTitle>
      <DialogContent>
        <Box sx={{ marginTop: 1 }}>
          {options.render({ values: values as TValues, setValues })}
        </Box>
      </DialogContent>
      <DialogActions>
        <Button type="button" onClick={(e) => {
          e.stopPropagation();
          close(false)
        }}>
          {options.noText}
        </Button>
        <Button type="button" variant="contained" onClick={(e) => {
          e.stopPropagation();
          close(true)
        }}>
          {options.yesText}
        </Button>
      </DialogActions>
    </Dialog>
  ) : null;

  return {
    openDialog,
    DialogComponent,
  } as const;
}

// ---------------- Example ----------------
// type Values = {
//   reason: string;
//   count: number;
// };
//
// const { openDialog, DialogComponent } = useDialog<Values>();
//
// const handleClick = async () => {
//   const result = await openDialog({
//     title: 'Confirm action',
//     initialValues: { reason: '', count: 1 },
//     render: ({ values, setValues }) => (
//       <>
//         <TextField
//           label="Reason"
//           fullWidth
//           value={values.reason}
//           onChange={(e) => setValues({...values, reason: e.target.value})}
//         />
//         <TextField
//           type="number"
//           label="Count"
//           value={values.count}
//           onChange={(e) => setValues({...values, count: Number(e.target.value)})}
//         />
//       </>
//     ),
//   });
//
//   if (result.confirmed) {
//     console.log(result.values);
//   }
// };
//
// return (
//   <>
//     <Button onClick={handleClick}>Open</Button>
//     {DialogComponent}
//   </>
// );
