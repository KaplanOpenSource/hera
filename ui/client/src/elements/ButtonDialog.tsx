import { ReactNode, useState } from 'react';
import { Dialog, DialogProps, IconButtonProps, Theme, ThemeProvider } from '@mui/material';
import { ButtonTooltip } from './ButtonTooltip';

export const ButtonDialog = ({
  icon,
  children,
  dialogProps,
  dialogTheme,
  onOpen,
  closeRef,
  ...props
}: Omit<IconButtonProps, 'children'> & {
  icon: ReactNode;
  title?: ReactNode;
  children: ReactNode | ((close: () => void) => ReactNode);
  dialogProps?: Omit<DialogProps, 'open' | 'onClose' | 'children'>;
  // Optional theme for the dialog only (trigger keeps the ambient theme). Use it to keep
  // the dialog readable when opened from a differently-themed area, e.g. the dark header.
  dialogTheme?: Theme;
  onOpen?: () => void;
  closeRef?: React.MutableRefObject<(() => void) | undefined>;
  // Render the trigger as a Button (icon + `title` label) instead of an IconButton.
  button?: boolean;
}) => {
  const [open, setOpen] = useState(false);

  const close = () => {
    setTimeout(() => {
      setOpen(false);
    }, 0);
  }

  if (closeRef) closeRef.current = close;

  const dialog = (
    <Dialog
      open={open}
      onClose={close}
      // The dialog is rendered inside other components (e.g. the tree); MUI portals it
      // to <body>, but React events still bubble up the React tree. Stop them so dialog
      // interactions don't reach the page underneath.
      onClick={(e) => e.stopPropagation()}
      onMouseDown={(e) => e.stopPropagation()}
      {...dialogProps}
    >
      {typeof children === 'function' ? children(close) : children}
    </Dialog>
  );

  return (
    <>
      <ButtonTooltip
        onClick={() => {
          onOpen?.();
          setOpen(true);
        }}
        {...props}
      >
        {icon}
      </ButtonTooltip>
      {dialogTheme ? <ThemeProvider theme={dialogTheme}>{dialog}</ThemeProvider> : dialog}
    </>
  );
};