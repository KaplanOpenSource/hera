import { Box, Button, ButtonProps, IconButton, IconButtonProps, Stack, Tooltip } from "@mui/material";
import { MouseEvent, ReactNode } from "react";

export const ButtonTooltip = ({
  onClick,
  children,
  title,
  button,
  ...restProps
}: {
  title?: ReactNode,
  onClick: NonNullable<IconButtonProps['onClick']>;
  children: any,
  // Render a regular Button (icon + `title` as label) instead of an IconButton.
  button?: boolean,
} & Omit<IconButtonProps, 'onClick'>) => {
  const handleClick = (e: MouseEvent<HTMLButtonElement>) => {
    e.stopPropagation();
    onClick(e);
  };

  return (
    <Tooltip title={title}>
      <Box>
        {button
          ? (
            <Button
              onClick={handleClick}
              size="small"
              startIcon={children}
              variant="contained"
              sx={{ borderRadius: 10, margin: '5px', alignItems: 'center' }}
              {...(restProps as ButtonProps)}
            >
              <span style={{marginBottom:'-5px'}}>
                {title}
              </span>
            </Button>
          )
          : (
            <IconButton
              onClick={handleClick}
              size="small"
              {...restProps}
            >
              {children}
            </IconButton>
          )}
      </Box>
    </Tooltip>
  );
}
