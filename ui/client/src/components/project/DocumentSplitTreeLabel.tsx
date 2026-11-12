import { Handyman, HotelClass } from "@mui/icons-material";
import { Typography, Stack, Tooltip } from "@mui/material";
import { SwitchCase, Case } from "../../elements/SwitchCase";
import { DESC_PATH_TOOLKIT, DESC_PATH_TYPE, VALUE_GROUP_REST, VALUE_GROUP_UNDEFINED } from "../../utils/splitTree";
import { DetailsViewToolkit } from "../details/toolkit/DetailsViewToolkit";
import { SplitTooltip } from "./SplitTooltip";
import { toolkitNameFromSplitId } from "../../shared/idDocId";
import { useProjectStore } from "../../stores/useProjectStore";
import { ProjectObj } from "../../objects/ProjectObj";

export const DocumentSplitTreeLabel = ({
  itemKey,
  path,
  value,
}: {
  itemKey: string,
  path: string,
  value: string,
}) => {
  // replacing '/item1/subitem2' to 'item1.subitem2'
  const fieldLabel = path.replace(/^\//, "").replace(/\//g, ".");
  const { currProject } = useProjectStore();
  const documents = currProject ? new ProjectObj(currProject).documents : [];
  const toolkitName = toolkitNameFromSplitId(itemKey, documents) ?? VALUE_GROUP_UNDEFINED;
  const isToolkit = path === DESC_PATH_TOOLKIT;


  return (
    <Tooltip
      title={isToolkit
        ? (
          <DetailsViewToolkit
            toolkitName={value === '' ? VALUE_GROUP_UNDEFINED : value}
          />
        ) : (
          <SplitTooltip
            fieldLabel={fieldLabel}
            value={value}
            toolkitName={toolkitName}
          />
        )
      }
      disableInteractive
      slotProps={{ tooltip: { sx: { maxWidth: 'none' } } }}
    >
      <span>
        <SwitchCase test={path}>
          <Case value={DESC_PATH_TOOLKIT}>
            <Stack direction={'row'} spacing={1}>
              <Handyman color="action" fontSize="small" />
              <b>
                {[VALUE_GROUP_UNDEFINED, ''].includes(value)
                  ? 'No toolkit'
                  : value
                }
              </b>
            </Stack>
          </Case>
          <Case value={DESC_PATH_TYPE}>
            <Stack direction={'row'} spacing={1}>
              <HotelClass color="action" fontSize="small" />
              <b>
                {[VALUE_GROUP_UNDEFINED, ''].includes(value)
                  ? 'No type'
                  : value
                }
              </b>
            </Stack>
          </Case>
          <Case isDefault={true}>
            <Typography>
              <SwitchCase test={value}>
                <Case value={VALUE_GROUP_REST}>
                  {fieldLabel} other values
                </Case>
                <Case value={VALUE_GROUP_UNDEFINED}>
                  {fieldLabel} not existing
                </Case>
                <Case isDefault={true}>
                  {fieldLabel} == <b>{value === '' ? '\'\'' : value}</b>
                </Case>
              </SwitchCase>
            </Typography>
          </Case>
        </SwitchCase>
      </span>
    </Tooltip>
  )
}
