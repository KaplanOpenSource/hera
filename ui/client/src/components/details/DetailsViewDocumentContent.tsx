import { Close, Done } from '@mui/icons-material';
import { Box, Divider, Stack, Typography } from '@mui/material';
import { SimpleTreeView } from '@mui/x-tree-view';
import { useEffect, useState } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { DocumentObj } from '../../objects/ProjectObj';
import { AgentConfig } from '../../shared/AgentConfig';
import { classifyDocument } from '../../shared/tabKind';
import { DocumentKindIcon } from '../project/DocumentKindIcon';
import { DATA_FORMAT_FIELD, DESC_FIELD, FORBIDDEN_FIELDS } from '../../shared/constants';
import { ProjectDocument, WorkflowDesc } from '../../shared/types';
import { getWorkflowSolver, isWorkflowDoc, setWorkflowSolver } from '../../shared/workflow';
import { copyOnly, copyWithout, reorderEntries } from '../../utils/utils';
import { AgentConfigEditor } from '../agents/AgentConfigEditor';
import { RunWorkflowButton } from '../workflow/RunWorkflowButton';
import { WorkflowEditor } from '../workflow/WorkflowEditor';
import { DeleteDocumentButton } from './DeleteDocumentButton';
import { DetailsViewDocumentHeader } from './DetailsViewDocumentHeader';
import { RawViewToggle } from './RawViewToggle';
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
  setDoc: (newDoc: DocumentObj) => Promise<void>,
  shownDoc: ProjectDocument,
  setShownDoc: (newDoc: ProjectDocument) => void,
}) => {
  const isAgent = isAgentConfigDoc(shownDoc);
  const isWorkflow = isWorkflowDoc(shownDoc);

  // OFF: the document's specialized view (agent/workflow editor, or the formulated
  // field tree). ON: the raw stored document with nothing hidden.
  const [rawView, setRawView] = useState(false);

  // Open each document in its specialized (non-raw) view.
  useEffect(() => {
    setRawView(false);
  }, [doc.docid]);

  const showFormulated = !rawView;
  const showAgentConfig = !rawView && isAgent;
  const showWorkflow = !rawView && isWorkflow;

  // The agent/workflow editors render the payload directly, so the unchangeable
  // meta fields move to the header (read-only) and are hidden from the editable
  // tree. Agent also hides resource (it IS the agent config); the workflow's
  // resource is a separate export path, so it stays editable.
  const showKindEditor = showAgentConfig || showWorkflow;
  let headerHiddenFields: string[] = [];
  if (showAgentConfig) {
    headerHiddenFields = ['resource', 'type', DATA_FORMAT_FIELD];
  } else if (showWorkflow) {
    headerHiddenFields = ['type', DATA_FORMAT_FIELD];
  }

  const isChanged = JSON.stringify(doc.data) !== JSON.stringify(shownDoc);
  return (
    <>
      <Stack direction={'row'} spacing={0.5} alignItems={'center'} justifyItems={'center'}>
        <DocumentKindIcon kind={classifyDocument(doc)} sx={{ color: 'text.secondary', mr: 1 }} />
        <Typography variant='h6' sx={{ marginRight: 1 }}>
          {doc.isConfig ? doc.project.name + ' config' : doc.name}
        </Typography>
        {isWorkflow && (
          <RunWorkflowButton
            projectName={doc.project.name}
            workflowName={(shownDoc.desc as WorkflowDesc).workflowName ?? doc.name}
            isChanged={isChanged}
            save={() => setDoc(new DocumentObj(shownDoc, doc.project))}
          />
        )}
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
        <Box sx={{ flexGrow: 1 }} />
        <RawViewToggle rawView={rawView} setRawView={setRawView} />
        <Divider orientation='vertical' flexItem sx={{ mx: 0.5 }} />
        <DeleteDocumentButton
          document={doc.data}
          projectName={doc.project.name}
          displayName={doc.isConfig ? doc.project.name + ' config' : doc.name}
          isConfig={doc.isConfig}
        />
      </Stack>
      <DetailsViewDocumentHeader
        docid={doc.docid}
        shownDoc={shownDoc}
        setShownDoc={setShownDoc}
        showFormulated={showFormulated}
        extraFields={!showKindEditor
          ? []
          : [
            { name: 'type', value: shownDoc.type },
            { name: DATA_FORMAT_FIELD, value: shownDoc.dataFormat },
          ]
        }
      />
      <SimpleTreeView
          defaultExpandedItems={[keyForDetailsViewItem(DESC_FIELD), keyForDetailsViewItem('resource')]}
          // Rows are near-symmetric now, so the chevron needs only a tiny nudge up
          // to line up with the name (required rows still reserve room below).
          sx={{ '& .MuiTreeItem-iconContainer': { transform: 'translateY(-1px)' } }}
        >
          {reorderEntries(Object.entries(shownDoc), [DESC_FIELD, 'resource']).map(([k, v]) => {
          if (FORBIDDEN_FIELDS.includes(k)) {
            return null;
          }
          if (headerHiddenFields.includes(k)) {
            return null;
          }
          const hideOnDesc = showFormulated && k === DESC_FIELD;
          const descHideFields = showWorkflow ? [...HIDE_ON_DESC, 'workflow', 'parameters'] : HIDE_ON_DESC;
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
                  // Rebuild from newVal (honors deletions) + only the hidden fields; merging over the full old desc would re-add deleted keys.
                  setShownDoc({ ...shownDoc, desc: { ...newVal, ...copyOnly(shownDoc.desc, descHideFields) } });
                }
              }}
            />
          );
        })}
        {showWorkflow && (
          <DetailsViewItem
            itemKey="solver"
            itemValue={getWorkflowSolver((shownDoc.desc as WorkflowDesc).workflow)}
            parentKey={undefined}
            setItemValue={newVal => setShownDoc({
              ...shownDoc,
              desc: { ...shownDoc.desc, workflow: setWorkflowSolver((shownDoc.desc as WorkflowDesc).workflow, newVal) } as WorkflowDesc,
            })}
          />
        )}
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
          actionButtons={
            <RunWorkflowButton
              projectName={doc.project.name}
              workflowName={(shownDoc.desc as WorkflowDesc).workflowName ?? doc.name}
              isChanged={isChanged}
              save={() => setDoc(new DocumentObj(shownDoc, doc.project))}
              sx={{ bgcolor: 'background.paper', boxShadow: 1, p: 0.5 }}
            />
          }
        />
      )}
    </>
  );
}
