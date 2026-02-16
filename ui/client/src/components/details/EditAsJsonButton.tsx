import { Keyboard } from "@mui/icons-material";
import { TextField } from "@mui/material";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { useDialog } from "../../elements/useDialog";

export const EditAsJsonButton = ({
  data,
  setData,
}: {
  data: any,
  setData: (newVal: any) => void,
}) => {
  const { DialogComponent, openDialog } = useDialog();
  return (
    <ButtonTooltip
      title={'Edit as Json Text'}
      onClick={async () => {
        const { confirmed, values } = await openDialog({
          title: 'Edit as Json Text',
          initialValues: { jsonText: JSON.stringify(data, null, 2), ok: true },
          render: ({ values, setValues }) => (
            <TextField
              label="Json"
              fullWidth
              multiline
              minRows={4}
              maxRows={20}
              value={values.jsonText}
              onChange={e => {
                e.stopPropagation();
                let ok = true;
                try { JSON.parse(e.target.value) } catch { ok = false }
                setValues({ ...values, jsonText: e.target.value, ok });
              }}
              error={!values.ok}
              helperText={values.ok ? '' : 'Invalid Json'}
            />
          )
        });
        if (confirmed && values && values.ok) {
          try {
            const json = JSON.parse(values.jsonText);
            setData(json);
          } catch { }
        }
      }}
    >
      <Keyboard />
      {DialogComponent}
    </ButtonTooltip>
  )
}