import { VisibilityOff } from '@mui/icons-material';
import { Box } from '@mui/material';
import { useEffect, useState } from 'react';
import { Group as PanelGroup, Panel, Separator as PanelResizeHandle } from 'react-resizable-panels';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { DocumentObj } from '../../objects/ProjectObj';
import { ProjectDocument } from '../../shared/types';
import { DetailsViewDocumentContent } from './DetailsViewDocumentContent';
import { PreviewChooser } from './PreviewChooser';

export const DetailsViewDocument = ({
  doc,
  setDoc,
  previewHidden = false,
  setPreviewHidden,
  previewAvailable = false,
}: {
  doc: DocumentObj,
  setDoc: (newDoc: DocumentObj) => void,
  previewHidden?: boolean,
  setPreviewHidden?: (hidden: boolean) => void,
  previewAvailable?: boolean,
}) => {
  const [shownDoc, setShownDoc] = useState<ProjectDocument>(JSON.parse(JSON.stringify(doc.data)));

  useEffect(() => {
    setShownDoc(JSON.parse(JSON.stringify(doc.data)));
  }, [doc.data])

  useEffect(() => {
    setPreviewHidden?.(false);
  }, [doc.docid])

  const showPreview = previewAvailable && !previewHidden;

  return showPreview
    ? (
      <Box sx={{ height: 'calc(100% + 32px)', m: -2 }}>
        <PanelGroup
          orientation="vertical"
          onLayoutChanged={(layout) => {
            if (layout['preview-panel'] === 0) {
              setPreviewHidden?.(true);
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
            id="preview-panel"
            defaultSize={50}
            minSize={5}
            collapsible
          >
            <Box sx={{ height: '100%', position: 'relative' }}>
              <Box sx={{ position: 'absolute', top: 4, right: 4, zIndex: 1000 }}>
                <ButtonTooltip
                  title="Hide preview"
                  onClick={() => setPreviewHidden?.(true)}
                  sx={{
                    backgroundColor: 'white',
                    '&:hover': { backgroundColor: '#eee' },
                  }}
                >
                  <VisibilityOff sx={{ fontSize: 14 }} />
                </ButtonTooltip>
              </Box>
              <PreviewChooser doc={shownDoc} />
            </Box>
          </Panel>
        </PanelGroup>
      </Box>
    )
    : (
      <DetailsViewDocumentContent
        doc={doc}
        setDoc={setDoc}
        shownDoc={shownDoc}
        setShownDoc={setShownDoc}
      />
    )
}
