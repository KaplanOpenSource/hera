import { Box } from '@mui/material';
import { ReactNode } from 'react';
import { Panel, Group as PanelGroup, Separator as PanelResizeHandle } from 'react-resizable-panels';

export const SplitWithSidebar = ({
  sidebar,
  collapsed,
  children,
}: {
  sidebar: ReactNode,
  collapsed: boolean,
  children: ReactNode,
}) => {
  return collapsed
    ? <Box sx={{ flex: 1, minWidth: 0 }}>{children}</Box>
    : (
      <PanelGroup orientation="horizontal">
        <Panel defaultSize={50} minSize={20}>
          {sidebar}
        </Panel>

        <PanelResizeHandle
          style={{
            width: 4,
            cursor: 'col-resize',
            backgroundColor: '#e0e0e0',
            outline: 'none',
          }}
        />

        <Panel defaultSize={50} minSize={20}>
          {children}
        </Panel>
      </PanelGroup>
    );
};
