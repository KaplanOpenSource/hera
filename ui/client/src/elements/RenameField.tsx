import { TextField } from "@mui/material";
import { useEffect, useState } from "react";
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

  useEffect(() => {
    setInternalValue(value);
  }, [value]);

  return (
    <TextField
      size='small'
      value={internalValue}
      onChange={(e) => setInternalValue(e.target.value)}
      onClick={e => e.stopPropagation()}
      onKeyDown={e => e.stopPropagation()}
      slotProps={{
        input: {
          endAdornment: (<>
            <ButtonTooltip
              title={'Rename'}
              onClick={() => setValue(internalValue)}
            >
              <Check />
            </ButtonTooltip>
            <ButtonTooltip
              title={'Cancel Rename'}
              onClick={() => setValue(value)}
            >
              <Close />
            </ButtonTooltip>
          </>)
        }
      }}
    />
  )
}