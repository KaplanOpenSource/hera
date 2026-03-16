import { Box } from '@mui/material';
import { useEffect, useState } from 'react';
import { Group as PanelGroup, Panel, Separator as PanelResizeHandle } from 'react-resizable-panels';
import { DocumentObj } from '../../objects/ProjectObj';
import { ProjectDocument } from '../../shared/types';
import { DetailsViewDocumentContent } from './DetailsViewDocumentContent';
import { isTileUrl, TileMapView } from './maps/TileMapView';

export const DetailsViewDocument = ({
  doc,
  setDoc,
}: {
  doc: DocumentObj,
  setDoc: (newDoc: DocumentObj) => void,
}) => {
  const [shownDoc, setShownDoc] = useState<ProjectDocument>(JSON.parse(JSON.stringify(doc.data)));

  useEffect(() => {
    setShownDoc(JSON.parse(JSON.stringify(doc.data)));
  }, [doc.data])

  const showMap = isTileUrl(shownDoc.resource);

  const detailsContent = (
    <DetailsViewDocumentContent
      doc={doc}
      setDoc={setDoc}
      shownDoc={shownDoc}
      setShownDoc={setShownDoc}
    />
  );

  if (!showMap) {
    return detailsContent;
  }

  return (
    <Box sx={{ height: 'calc(100% + 32px)', m: -2 }}>
      <PanelGroup orientation="vertical">
        <Panel defaultSize={50} minSize={20}>
          <Box sx={{ height: '100%', overflow: 'auto', p: 2 }}>
            {detailsContent}
          </Box>
        </Panel>

        <PanelResizeHandle
          style={{
            height: 4,
            cursor: 'row-resize',
            backgroundColor: '#e0e0e0',
            outline: 'none',
          }}
        />

        <Panel defaultSize={50} minSize={20}>
          <TileMapView url={shownDoc.resource as string} />
        </Panel>
      </PanelGroup>
    </Box>
  )
}
