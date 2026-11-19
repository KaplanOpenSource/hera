import { FormControlLabel, Switch } from '@mui/material';

// A TRUE/FALSE switch for a boolean field.
export const DetailsViewItemBoolean = ({
  itemValue,
  setItemValue,
}: {
  itemValue: boolean,
  setItemValue: (newVal: boolean) => void,
}) => {
  let label = 'FALSE';
  if (itemValue) {
    label = 'TRUE';
  }

  return (
    <FormControlLabel
      label={label}
      onClick={e => e.stopPropagation()}
      onKeyDown={e => e.stopPropagation()}
      control={
        <Switch
          size='small'
          checked={itemValue}
          onChange={(e) => setItemValue(e.target.checked)}
        />
      }
    />
  );
};
