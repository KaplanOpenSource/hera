import { Box, IconButton, IconButtonProps, Tooltip } from "@mui/material";
import { MouseEvent, ReactNode } from "react";

export const ButtonTooltip = ({
  onClick,
  children,
  title,
  ...restProps
}: {
  title?: ReactNode,
  onClick: NonNullable<IconButtonProps['onClick']>;
  children: any,
} & Omit<IconButtonProps, 'onClick'>) => {
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