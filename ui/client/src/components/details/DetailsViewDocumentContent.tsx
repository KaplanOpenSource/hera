import { AccountTree, Close, Done, DynamicForm, Science } from '@mui/icons-material';
import { Stack, Typography } from '@mui/material';
import { SimpleTreeView } from '@mui/x-tree-view';
import { useEffect, useState } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { DocumentObj } from '../../objects/ProjectObj';
import { AgentConfig } from '../../shared/AgentConfig';
import { FORBIDDEN_FIELDS } from '../../shared/constants';
import { ProjectDocument, WorkflowDesc } from '../../shared/types';
import { isWorkflowDoc } from '../../shared/workflow';
import { copyWithout, reorderEntries } from '../../utils/utils';
import { AgentConfigEditor } from '../agents/AgentConfigEditor';
import { WorkflowEditor } from '../workflow/WorkflowEditor';
import { DetailsViewDocumentHeader } from './DetailsViewDocumentHeader';
import { DetailsViewItem, keyForDetailsViewItem } from './DetailsViewItem';

const HIDE_ON_DESC = ['datasourceName', 'toolkit', 'version'];
const isAgentConfigDoc = (doc: ProjectDocument) => {
  return doc && typeof doc?.resource === 'object' && doc?.resource.effects !== undefined;
}

export const DetailsViewDocumentContent = ({
  doc,
  setDoc,
  shownDoc,
  setShownDoc,
}: {
  doc: DocumentObj,
  setDoc: (newDoc: DocumentObj) => void,
  shownDoc: ProjectDocument,
  setShownDoc: (newDoc: ProjectDocument) => void,
}) => {
  const [showFormulated, setShowFormulated] = useState(true);
  const [showAgentConfig, setShowAgentConfig] = useState(isAgentConfigDoc(doc.data));
  const [showWorkflow, setShowWorkflow] = useState(isWorkflowDoc(doc.data));

  const isAgent = isAgentConfigDoc(shownDoc);
  const isWorkflow = isWorkflowDoc(shownDoc);

  useEffect(() => {
    if (!isAgent) {
      setShowAgentConfig(false);
    }
  }, [isAgent]);

  useEffect(() => {
    if (!isWorkflow) {
      setShowWorkflow(false);
    }
  }, [isWorkflow]);

  useEffect(() => {
    setShowAgentConfig(isAgentConfigDoc(doc.data));
    setShowWorkflow(isWorkflowDoc(doc.data));
  }, [doc.docid]);

  useEffect(() => {
    if (!showFormulated) {
      setShowAgentConfig(false);
      setShowWorkflow(false);
    }
  }, [showFormulated])

  useEffect(() => {
    if (showAgentConfig) {
      setShowFormulated(true);
    }
  }, [showAgentConfig])

  useEffect(() => {
    if (showWorkflow) {
      setShowFormulated(true);
    }
  }, [showWorkflow])

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
        <ButtonTooltip
          title={'Show Workflow'}
          onClick={() => setShowWorkflow(!showWorkflow)}
          disabled={!isWorkflow}
        >
          <AccountTree color={showWorkflow ? 'primary' : 'inherit'} />
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
          const descHideFields = showWorkflow ? [...HIDE_ON_DESC, 'workflow'] : HIDE_ON_DESC;
          return (
            <DetailsViewItem
              key={k}
              itemKey={k}
              itemValue={!hideOnDesc ? v : copyWithout(v, descHideFields)}
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
      {showWorkflow && (
        <WorkflowEditor
          workflow={(shownDoc.desc as WorkflowDesc).workflow}
          setWorkflow={newVal => setShownDoc({ ...shownDoc, desc: { ...shownDoc.desc, workflow: newVal } as WorkflowDesc })}
        />
      )}
    </>
  );
}
