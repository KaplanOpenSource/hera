import { FormControl, InputLabel, Select, SelectChangeEvent, MenuItem } from "@mui/material"

export const SelectProperty = ({
  label,
  value,
  setValue,
  menuItems,
}: {
  label: string,
  value: string,
  setValue: (v: string) => void,
  menuItems: { name: string }[],
}) => {
  return (
    <FormControl style={{ minWidth: '110px' }}>
      <InputLabel>
        {label}
      </InputLabel>
      <Select
        value={value}
        label={label}
        size='small'
        onChange={(e: SelectChangeEvent) => setValue(e.target.value)}
      >
        {menuItems.map(({ name }) => (
          <MenuItem key={name}
            value={name}
          >
            {name}
          </MenuItem>
        ))}
      </Select>
    </FormControl>
  )
}