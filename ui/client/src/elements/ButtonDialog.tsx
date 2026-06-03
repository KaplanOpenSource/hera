import { ReactNode, useState } from 'react';
import { Dialog, DialogProps, IconButtonProps } from '@mui/material';
import { ButtonTooltip } from './ButtonTooltip';

export const ButtonDialog = ({
  icon,
  children,
  dialogProps,
  onOpen,
  closeRef,
  ...props
}: Omit<IconButtonProps, 'children'> & {
  icon: ReactNode;
  title?: ReactNode;
  children: ReactNode | ((close: () => void) => ReactNode);
  dialogProps?: Omit<DialogProps, 'open' | 'onClose' | 'children'>;
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
      <Dialog
        open={open}
        onClose={close}
        {...dialogProps}
      >
        {typeof children === 'function' ? children(close) : children}
      </Dialog>
    </>
  );
};