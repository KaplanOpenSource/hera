import { FormGroup, FormControlLabel, Checkbox } from "@mui/material"

export const BooleanProperty = ({
  label,
  value,
  setValue,
}: {
  label: string,
  value: boolean,
  setValue: (val: boolean) => void,
}) => {
  return (
    <FormGroup>
      <FormControlLabel
        label={label}
        control={<Checkbox
          checked={value}
          onChange={(e) => {
            e.stopPropagation();
            setValue(e.target.checked);
          }}
        />}
      />
    </FormGroup>
  )
}