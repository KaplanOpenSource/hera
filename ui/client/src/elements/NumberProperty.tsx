import { TextField, TextFieldProps } from "@mui/material";

export const NumberProperty = ({
  value,
  setValue,
  min,
  max,
  step,
  ...props
}: Omit<TextFieldProps, 'size' | 'variant' | 'value' | 'onChange' | 'type'> & {
  value: number;
  setValue: (v: number) => void;
  min?: number;
  max?: number;
  step?: number;
}) => (
  <TextField
    size="small"
    variant="outlined"
    type="number"
    value={value}
    onClick={(e) => e.stopPropagation()}
    onChange={(e) => {
      e.stopPropagation();
      setValue(Number(e.target.value));
    }}
    slotProps={{
      htmlInput: { min, max, step },
    }}
    {...props}
  />
);