import { Check } from "@mui/icons-material";
import { Box, TextField } from "@mui/material";
import { ButtonTooltip } from "../../elements/ButtonTooltip";

export const TextFieldWithApply = ({
  value,
  setValue,
  onApply,
  textTitle,
  buttonTitle,
}: {
  value: string,
  setValue: (v: string) => void,
  onApply: () => void,
  textTitle: string,
  buttonTitle: string,
}) => {
  return (
    <Box position={'relative'}>
      <TextField
        fullWidth
        label={textTitle}
        value={value}
        onClick={e => e.stopPropagation()}
        onChange={(e) => {
          setValue(e.target.value);
        }}
        size="small"
        rows={15}
        multiline={true}
      />
      <ButtonTooltip
        sx={{
          position: 'absolute',
          top: '8px',
          right: '20px',
          zIndex: 1
        }}
        title={buttonTitle}
        onClick={() => onApply()}
      >
        <Check />
      </ButtonTooltip>
    </Box>
  )
}