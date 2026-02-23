import { ReactNode, useState } from 'react';
import { Dialog, DialogProps, IconButtonProps } from '@mui/material';
import { ButtonTooltip } from './ButtonTooltip';

export const ButtonDialog = ({
  icon,
  children,
  dialogProps,
  ...props
}: Omit<IconButtonProps, 'children'> & {
  icon: ReactNode;
  title?: ReactNode;
  children: ReactNode;
  dialogProps?: Omit<DialogProps, 'open' | 'onClose' | 'children'>;
}) => {
  const [open, setOpen] = useState(false);

  return (
    <>
      <ButtonTooltip
        onClick={() => setOpen(true)}
        {...props}
      >
        {icon}
      </ButtonTooltip>
      <Dialog
        open={open}
        onClose={() => setOpen(false)}
        {...dialogProps}
      >
        {children}
      </Dialog>
    </>
  );
};