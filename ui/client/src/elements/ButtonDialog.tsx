import { ReactNode, useState } from 'react';
import { Dialog, DialogProps, IconButtonProps } from '@mui/material';
import { ButtonTooltip } from './ButtonTooltip';

export const ButtonDialog = ({
  icon,
  children,
  dialogProps,
  onOpen,
  ...props
}: Omit<IconButtonProps, 'children'> & {
  icon: ReactNode;
  title?: ReactNode;
  children: ReactNode | ((close: () => void) => ReactNode);
  dialogProps?: Omit<DialogProps, 'open' | 'onClose' | 'children'>;
  onOpen?: () => void;
}) => {
  const [open, setOpen] = useState(false);
  const close = () => setOpen(false);

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