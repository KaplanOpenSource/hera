import { Handyman } from "@mui/icons-material";
import { Typography, Stack } from "@mui/material";
import { SwitchCase, Case } from "../../elements/SwitchCase";
import { DESC_PATH_TOOLKIT, VALUE_GROUP_REST, VALUE_GROUP_UNDEFINED } from "./DocumentSplitGroup";

export const DocumentSplitTreeLabel = ({
  path,
  value,
}: {
  path: string,
  value: string,
}) => {
  // replacing '/item1/subitem2' to 'item1.subitem2'
  const fieldLabel = path.replace(/^\//, "").replace(/\//g, ".");

  return (
    <Typography>
      <SwitchCase test={path}>
        <Case value={DESC_PATH_TOOLKIT}>
          <Stack direction={'row'} spacing={1}>
            <Handyman />
            <b>
              {[VALUE_GROUP_UNDEFINED].includes(value)
                ? 'No toolkit'
                : value
              }
            </b>
          </Stack>
        </Case>
        <Case isDefault={true}>
          <SwitchCase test={value}>
            <Case value={VALUE_GROUP_REST}>
              <b>{fieldLabel}</b> other values
            </Case>
            <Case value={VALUE_GROUP_UNDEFINED}>
              <b>{fieldLabel}</b> not existing
            </Case>
            <Case isDefault={true}>
              <b>{fieldLabel}</b> == <b>{value === '' ? '\'\'' : value}</b>
            </Case>
          </SwitchCase>
        </Case>
      </SwitchCase>
    </Typography>
  )
}
