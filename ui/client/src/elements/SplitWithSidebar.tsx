import { ReactNode, useEffect } from 'react';
import { Panel, Group as PanelGroup, Separator as PanelResizeHandle, usePanelRef } from 'react-resizable-panels';

export const SplitWithSidebar = ({
  sidebar,
  collapsed,
  children,
}: {
  sidebar: ReactNode,
  collapsed: boolean,
  children: ReactNode,
}) => {
  const sidebarRef = usePanelRef();

  useEffect(() => {
    if (collapsed) {
      sidebarRef.current?.collapse();
    } else {
      sidebarRef.current?.expand();
    }
  }, [collapsed]);

  return (
    <PanelGroup orientation="horizontal">
      <Panel
        panelRef={sidebarRef}
        defaultSize={50}
        minSize={20}
        collapsible
      >
        {sidebar}
      </Panel>

      <PanelResizeHandle
        style={{
          width: collapsed ? 0 : 4,
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
