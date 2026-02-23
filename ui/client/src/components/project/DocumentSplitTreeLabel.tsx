import { Handyman } from "@mui/icons-material";
import { Typography, Stack } from "@mui/material";
import { SwitchCase, Case } from "../../elements/SwitchCase";
import { TOOLKIT_DESC_PATH, VALUE_GROUP_REST, VALUE_GROUP_UNDEFINED } from "./DocumentSplitGroup";

export const DocumentSplitTreeLabel = ({
  path,
  value,
}: {
  path: string,
  value: string,
}) => {
  const isToolkit = path === TOOLKIT_DESC_PATH;

  // replacing '/item1/subitem2' to 'item1.subitem2'
  const fieldLabel = path.replace(/^\//, "").replace(/\//g, ".");

  return (
    <Typography>
      {isToolkit
        ? (
          <Stack direction={'row'} spacing={1}>
            <Handyman />
            <b>
              {value === VALUE_GROUP_UNDEFINED ? 'No toolkit' : value}
            </b>
          </Stack>
        )
        : (
          <SwitchCase test={value}>
            <Case value={VALUE_GROUP_REST}>
              <b>{fieldLabel}</b> other values
            </Case>
            <Case value={VALUE_GROUP_UNDEFINED}>
              <b>{fieldLabel}</b> not existing
            </Case>
            <Case isDefault={true}>
              <b>{fieldLabel}</b> == <b>{value}</b>
            </Case>
          </SwitchCase>
        )}
    </Typography>
  )
}
