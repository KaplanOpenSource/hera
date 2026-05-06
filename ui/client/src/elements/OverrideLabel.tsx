import { Tooltip, Typography } from "@mui/material";
import { ReactNode } from "react";
import { OverrideTooltipTable } from "./OverrideTooltipTable";

export const OverrideLabel = ({
  children,
  isOverridden,
  treePathKey,
  overrides,
  repoJsons,
}: {
  children: ReactNode;
  isOverridden: boolean;
  treePathKey: string;
  overrides?: { [path: string]: string[] };
  repoJsons?: { [path: string]: { [key: string]: any } };
}) => {
  if (isOverridden && repoJsons && overrides) {
    return (
      <Tooltip
        slotProps={{ tooltip: { sx: { maxWidth: 'none' } } }}
        title={
          <OverrideTooltipTable
            overridePath={treePathKey}
            filePaths={overrides[treePathKey]}
            repoJsons={repoJsons}
          />
        }
      >
        <Typography variant="body2" color="error">{children}</Typography>
      </Tooltip>
    );
  }

  return <Typography variant="body2">{children}</Typography>;
};
