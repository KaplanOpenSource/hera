import { Map } from '@mui/icons-material';
import { Box } from '@mui/material';
import { useEffect, useState } from 'react';
import { Group as PanelGroup, Panel, Separator as PanelResizeHandle } from 'react-resizable-panels';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
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
  const [mapHidden, setMapHidden] = useState(false);

  useEffect(() => {
    setShownDoc(JSON.parse(JSON.stringify(doc.data)));
  }, [doc.data])

  const hasMap = isTileUrl(shownDoc.resource);
  const showMap = hasMap && !mapHidden;

  return showMap
    ? (
      <Box sx={{ height: 'calc(100% + 32px)', m: -2 }}>
        <PanelGroup
          orientation="vertical"
          onLayoutChanged={(layout) => {
            if (layout['map-panel'] === 0) {
              setMapHidden(true);
            }
          }}
        >
          <Panel defaultSize={50} minSize={20}>
            <Box sx={{ height: '100%', overflow: 'auto', p: 2 }}>
              <DetailsViewDocumentContent
                doc={doc}
                setDoc={setDoc}
                shownDoc={shownDoc}
                setShownDoc={setShownDoc}
              />
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

          <Panel
            id="map-panel"
            defaultSize={50}
            minSize={5}
            collapsible
          >
            <TileMapView
              url={shownDoc.resource as string}
              onClose={() => setMapHidden(true)}
            />
          </Panel>
        </PanelGroup>
      </Box>
    )
    : (
      <Box sx={{ height: '100%', position: 'relative' }}>
        <DetailsViewDocumentContent
          doc={doc}
          setDoc={setDoc}
          shownDoc={shownDoc}
          setShownDoc={setShownDoc}
        />
        {hasMap && mapHidden && (
          <Box sx={{ position: 'absolute', bottom: 4, right: 4 }}>
            <ButtonTooltip
              title="Show map"
              onClick={() => setMapHidden(false)}
            >
              <Map sx={{ fontSize: 14 }} />
            </ButtonTooltip>
          </Box>
        )}
      </Box>
    )
}
