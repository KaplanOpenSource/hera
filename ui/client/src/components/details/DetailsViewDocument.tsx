import { Close, Done, DynamicForm, Science } from '@mui/icons-material';
import { Stack, Typography } from '@mui/material';
import { SimpleTreeView } from '@mui/x-tree-view';
import { useEffect, useState } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { DocumentObj } from '../../objects/ProjectObj';
import { AgentConfig } from '../../shared/AgentConfig';
import { FORBIDDEN_FIELDS } from '../../shared/constants';
import { ProjectDocument } from '../../shared/types';
import { copyWithout, reorderEntries } from '../../utils/utils';
import { AgentConfigEditor } from '../agents/AgentConfigEditor';
import { DetailsViewDocumentHeader } from './DetailsViewDocumentHeader';
import { DetailsViewItem, keyForDetailsViewItem } from './DetailsViewItem';
import { isTileUrl, TileMapView } from './TileMapView';

const HIDE_ON_DESC = ['datasourceName', 'toolkit', 'version'];
const isAgentConfigDoc = (doc: ProjectDocument) => {
  return doc && typeof doc?.resource === 'object' && doc?.resource.effects !== undefined;
}

export const DetailsViewDocument = ({
  doc,
  setDoc,
}: {
  doc: DocumentObj,
  setDoc: (newDoc: DocumentObj) => void,
}) => {
  const [shownDoc, setShownDoc] = useState<ProjectDocument>(JSON.parse(JSON.stringify(doc.data)));
  const [showFormulated, setShowFormulated] = useState(true);
  const [showAgentConfig, setShowAgentConfig] = useState(isAgentConfigDoc(doc.data));

  const isAgent = isAgentConfigDoc(shownDoc);

  useEffect(() => {
    if (!isAgent) {
      setShowAgentConfig(false);
    }
  }, [isAgent]);

  useEffect(() => {
    setShowAgentConfig(isAgentConfigDoc(doc.data));
  }, [doc.docid]);

  useEffect(() => {
    setShownDoc(JSON.parse(JSON.stringify(doc.data)));
  }, [doc.data])

  useEffect(() => {
    if (!showFormulated) {
      setShowAgentConfig(false);
    }
  }, [showFormulated])

  useEffect(() => {
    if (showAgentConfig) {
      setShowFormulated(true);
    }
  }, [showAgentConfig])

  const isChanged = JSON.stringify(doc.data) !== JSON.stringify(shownDoc);
  return (
    <>
      <Stack direction={'row'} alignItems={'center'} justifyItems={'center'}>
        <Typography variant='h6' sx={{ marginRight: 1 }}>
          {doc.isConfig ? doc.project.name + ' config' : doc.name}
        </Typography>
        <ButtonTooltip
          title={'Show Formulated'}
          onClick={() => setShowFormulated(!showFormulated)}
        >
          <DynamicForm color={showFormulated ? 'primary' : 'inherit'} />
        </ButtonTooltip>
        <ButtonTooltip
          title={'Show Agent Config'}
          onClick={() => setShowAgentConfig(!showAgentConfig)}
          disabled={!isAgent}
        >
          <Science color={showAgentConfig ? 'primary' : 'inherit'} />
        </ButtonTooltip>
        {isChanged
          ? (<>
            <ButtonTooltip
              title='Update Document'
              onClick={() => setDoc(new DocumentObj(shownDoc, doc.project))}
            >
              <Done />
            </ButtonTooltip>
            <ButtonTooltip
              title='Revert Document'
              onClick={() => setShownDoc(JSON.parse(JSON.stringify(doc.data)))}
            >
              <Close />
            </ButtonTooltip>
          </>)
          : null}
      </Stack>
      <DetailsViewDocumentHeader
        docid={doc.docid}
        shownDoc={shownDoc}
        setShownDoc={setShownDoc}
        showFormulated={showFormulated}
        extraFields={!showAgentConfig
          ? []
          : [
            { name: 'type', value: shownDoc.type },
            { name: 'dataFormat', value: shownDoc.dataFormat },
          ]
        }
      />
      <SimpleTreeView
        defaultExpandedItems={[keyForDetailsViewItem('desc'), keyForDetailsViewItem('resource')]}
      >
        {reorderEntries(Object.entries(shownDoc), ['desc', 'resource']).map(([k, v]) => {
          if (FORBIDDEN_FIELDS.includes(k)) {
            return null;
          }
          if (showAgentConfig && ['resource', 'type', 'dataFormat'].includes(k)) {
            return null;
          }
          const hideOnDesc = showFormulated && k === 'desc';
          return (
            <DetailsViewItem
              key={k}
              itemKey={k}
              itemValue={!hideOnDesc ? v : copyWithout(v, HIDE_ON_DESC)}
              parentKey={undefined}
              setItemValue={newVal => {
                if (!hideOnDesc) {
                  setShownDoc({ ...shownDoc, [k]: newVal });
                } else {
                  setShownDoc({ ...shownDoc, desc: { ...shownDoc.desc, ...newVal } });
                }
              }}
            />
          );
        })}
      </SimpleTreeView>
      {showAgentConfig
        ? (
          <AgentConfigEditor
            agentResource={shownDoc.resource as AgentConfig}
            setAgentResource={newVal => setShownDoc({ ...shownDoc, resource: newVal })}
          />
        )
        : null
      }
      {isTileUrl(shownDoc.resource)
        ? <TileMapView url={shownDoc.resource} />
        : null
      }
    </>
  )
}
