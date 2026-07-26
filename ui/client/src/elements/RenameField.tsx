import { TextField } from "@mui/material";
import { useEffect, useRef, useState } from "react";
import { ButtonTooltip } from "./ButtonTooltip";
import { Check, Close } from "@mui/icons-material";

export const RenameField = ({
  value,
  setValue,
}: {
  value: string,
  setValue: (newVal: string) => void,
}) => {
  const [internalValue, setInternalValue] = useState(value);
  const inputRef = useRef<HTMLInputElement>(null);

  useEffect(() => {
    setInternalValue(value);
  }, [value]);

  useEffect(() => {
    const input = inputRef.current;
    if (input) {
      input.focus();
      input.setSelectionRange(input.value.length, input.value.length);
    }
  }, []);

  return (
    <TextField
      size='small'
      inputRef={inputRef}
      value={internalValue}
      onChange={(e) => setInternalValue(e.target.value)}
      onClick={e => e.stopPropagation()}
      onKeyDown={e => {
        e.stopPropagation();
        if (e.key === 'Enter') {
          setValue(internalValue);
        } else if (e.key === 'Escape') {
          setValue(value);
        }
      }}
      slotProps={{
        input: {
          endAdornment: (<>
            <ButtonTooltip
              title={'Rename'}
              onClick={() => setValue(internalValue)}
            >
              <Check fontSize="small" />
            </ButtonTooltip>
            <ButtonTooltip
              title={'Cancel Rename'}
              onClick={() => setValue(value)}
            >
              <Close fontSize="small" />
            </ButtonTooltip>
          </>)
        }
      }}
      fullWidth
    />
  )
}