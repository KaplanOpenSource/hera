import { TextField, Typography } from "@mui/material";
import { useEffect, useRef, useState } from "react";
import { ButtonTooltip } from "./ButtonTooltip";
import { Check, Close } from "@mui/icons-material";

export const RenameField = ({
  value,
  setValue = undefined,
  defaultEditing = false,
  labelMinWidth = undefined,
}: {
  value: string,
  setValue?: (newVal: string) => void,
  defaultEditing?: boolean,
  labelMinWidth?: string,
}) => {
  const [editing, setEditing] = useState(defaultEditing);
  const [internalValue, setInternalValue] = useState(value);
  const inputRef = useRef<HTMLInputElement>(null);

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

  const confirm = () => {
    setValue?.(internalValue);
    setEditing(false);
  };

  const cancel = () => {
    setInternalValue(value);
    setEditing(false);
  };

  return (
    (editing && setValue)
      ? (
        <TextField
          size='small'
          inputRef={inputRef}
          value={internalValue}
          onChange={(e) => setInternalValue(e.target.value)}
          onClick={e => e.stopPropagation()}
          onKeyDown={e => {
            e.stopPropagation();
            if (e.key === 'Enter') {
              confirm();
            } else if (e.key === 'Escape') {
              cancel();
            }
          }}
          slotProps={{
            input: {
              endAdornment: (<>
                <ButtonTooltip title={'Rename'} onClick={confirm}>
                  <Check fontSize="small" />
                </ButtonTooltip>
                <ButtonTooltip title={'Cancel Rename'} onClick={cancel}>
                  <Close fontSize="small" />
                </ButtonTooltip>
              </>)
            }
          }}
        />
      )
      : (
        <Typography
          onClick={() => setValue && setEditing(true)}
          sx={{
            whiteSpace: 'nowrap',
            minWidth: labelMinWidth,
            cursor: setValue ? 'text' : 'default'
          }}
        >
          {value}
        </Typography>
      )
  )
}
