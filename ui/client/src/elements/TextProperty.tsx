import { TextField, TextFieldProps } from '@mui/material';

export const TextProperty = ({
  value,
  setValue,
  ...props
}: Omit<TextFieldProps, 'size' | 'variant' | 'value' | 'onChange'> & {
  value: string;
  setValue: (v: string) => void;
}) => (
  <TextField
    size="small"
    variant="outlined"
    value={value}
    onClick={(e) => e.stopPropagation()}
    onChange={(e) => {
      e.stopPropagation();
      setValue(e.target.value);
    }}
    {...props}
  />
);