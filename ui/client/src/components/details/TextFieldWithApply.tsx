import { Check, UploadFile } from "@mui/icons-material";
import { Box, Stack, TextField } from "@mui/material";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { useRef } from "react";

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
  const inputRef = useRef<HTMLInputElement | null>(null);

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
      <Stack
        direction={'row'}
        sx={{
          position: 'absolute',
          top: '8px',
          right: '20px',
          zIndex: 1
        }}
      >
        <ButtonTooltip
          title={'Load text file'}
          onClick={() => inputRef.current?.click()}
        >
          <UploadFile />
          <input
            ref={inputRef}
            type="file"
            style={{ display: 'none' }}
            accept=".txt,.json"
            onChange={async (e) => {
              const file = (e.target.files || [])[0];
              if (file) {
                const str = await file.text();
                setValue(str);
              }
            }}
          />
        </ButtonTooltip>
        <ButtonTooltip
          title={buttonTitle}
          onClick={() => onApply()}
        >
          <Check />
        </ButtonTooltip>
      </Stack>
    </Box>
  )
}