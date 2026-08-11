import { Box, TextField, Typography } from "@mui/material";
import { ReactNode, useEffect, useRef, useState } from "react";

export const RenameField = ({
  value,
  setValue = undefined,
  defaultEditing = false,
  labelMinWidth = undefined,
  valueForView = undefined,
}: {
  value: string,
  setValue?: (newVal: string) => void,
  defaultEditing?: boolean,
  labelMinWidth?: string,
  // What to render in view mode instead of `value` (editing still uses `value`).
  valueForView?: ReactNode,
}) => {
  const [editing, setEditing] = useState(defaultEditing);
  const [internalValue, setInternalValue] = useState(value);
  const inputRef = useRef<HTMLInputElement>(null);
  const cancelOnBlur = useRef(false);

  useEffect(() => {
    setInternalValue(value);
  }, [value]);

  useEffect(() => {
    if (editing) {
      const input = inputRef.current;
      if (input) {
        input.focus();
        input.setSelectionRange(input.value.length, input.value.length);
      }
    }
  }, [editing]);

  // Blur is the single save path: Enter blurs to save, Escape sets cancelOnBlur then blurs.
  const handleBlur = () => {
    if (cancelOnBlur.current) {
      cancelOnBlur.current = false;
      setInternalValue(value);
    } else {
      setValue?.(internalValue);
    }
    setEditing(false);
  };

  const viewLabel = valueForView !== undefined
    ? (
      <Box
        onClick={() => setValue && setEditing(true)}
        sx={{
          display: 'flex',
          minWidth: 0,
          overflow: 'hidden',
          cursor: setValue ? 'text' : 'default',
        }}
      >
        {valueForView}
      </Box>
    )
    : (
      <Typography
        onClick={() => setValue && setEditing(true)}
        sx={{
          whiteSpace: 'nowrap',
          minWidth: labelMinWidth,
          flexShrink: 0,
          cursor: setValue ? 'text' : 'default'
        }}
      >
        {value}
      </Typography>
    );

  return (
    (editing && setValue)
      ? (
        <TextField
          size='small'
          inputRef={inputRef}
          value={internalValue}
          onChange={(e) => setInternalValue(e.target.value)}
          onClick={e => e.stopPropagation()}
          onBlur={handleBlur}
          onKeyDown={e => {
            e.stopPropagation();
            if (e.key === 'Enter') {
              inputRef.current?.blur();
            } else if (e.key === 'Escape') {
              cancelOnBlur.current = true;
              inputRef.current?.blur();
            }
          }}
        />
      )
      : viewLabel
  )
}
