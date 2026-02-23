import { Box, IconButton, IconButtonProps, Tooltip } from "@mui/material";
import { MouseEvent, ReactNode } from "react";

export const ButtonTooltip = ({
  onClick,
  children,
  title,
  ...restProps
}: {
  onClick: (e: MouseEvent<HTMLElement>) => void,
  title?: ReactNode,
  children: any,
} & IconButtonProps) => {
  return (
    <Tooltip
      title={title}
    >
      <Box>
        <IconButton
          onClick={e => {
            e.stopPropagation();
            onClick(e);
          }}
          size="small"
          {...restProps}
        >
          {children}
        </IconButton>
      </Box>
    </Tooltip>
  )
}