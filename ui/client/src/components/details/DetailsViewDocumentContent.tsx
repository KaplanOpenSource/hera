import { Close, Done } from '@mui/icons-material';
import { Stack, Typography } from '@mui/material';
import { SimpleTreeView } from '@mui/x-tree-view';
import { useEffect, useState } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { DocumentObj } from '../../objects/ProjectObj';
import { AgentConfig } from '../../shared/AgentConfig';
import { DATA_FORMAT_FIELD, DESC_FIELD, FORBIDDEN_FIELDS } from '../../shared/constants';
import { TabKind } from '../../shared/tabKind';
import { ProjectDocument, WorkflowDesc } from '../../shared/types';
import { getWorkflowSolver, isWorkflowDoc, setWorkflowSolver } from '../../shared/workflow';
import { copyOnly, copyWithout, reorderEntries } from '../../utils/utils';
import { AgentConfigEditor } from '../agents/AgentConfigEditor';
import { RunWorkflowButton } from '../workflow/RunWorkflowButton';
import { WorkflowEditor } from '../workflow/WorkflowEditor';
import { DetailsViewDocumentHeader } from './DetailsViewDocumentHeader';
import { DocView, DocViewSelector } from './DocViewSelector';
import { DetailsViewItem, keyForDetailsViewItem } from './DetailsViewItem';
import { DetailsVisibility, DetailsVisibilityToggle } from './DetailsVisibilityToggle';

const HIDE_ON_DESC = ['datasourceName', 'toolkit', 'version'];
const isAgentConfigDoc = (doc: ProjectDocument) => {
  return doc && typeof doc?.resource === 'object' && doc?.resource.effects !== undefined;
}

const defaultView = (doc: ProjectDocument): DocView => {
  return isAgentConfigDoc(doc) ? TabKind.Agent : isWorkflowDoc(doc) ? TabKind.Workflow : 'formulated';
};

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
  const isAgent = isAgentConfigDoc(shownDoc);
  const isWorkflow = isWorkflowDoc(shownDoc);

  const [docView, setDocView] = useState<DocView>(() => defaultView(doc.data));
  const [detailsVisibility, setDetailsVisibility] = useState<DetailsVisibility>(DetailsVisibility.Both);
  const showHeader = detailsVisibility === DetailsVisibility.Both;
  const showTree = detailsVisibility !== DetailsVisibility.None;

  // When switching to a different document, reset to its default view.
  useEffect(() => {
    setDocView(defaultView(doc.data));
  }, [doc.docid]);

  // The agent/workflow views only apply to those document kinds; if the shown
  // doc is no longer one of them, fall back to the formulated view.
  useEffect(() => {
    if ((docView === TabKind.Agent && !isAgent) || (docView === TabKind.Workflow && !isWorkflow)) {
      setDocView('formulated');
    }
  }, [isAgent, isWorkflow, docView]);

  const showFormulated = docView !== 'raw';
  const showAgentConfig = docView === TabKind.Agent;
  const showWorkflow = docView === TabKind.Workflow;

  // Agent/workflow views render the document's payload in a dedicated editor, so
  // the unchangeable meta fields move to the header (read-only) and are hidden
  // from the editable tree. Agent also hides resource (it IS the agent config);
  // the workflow's resource is a separate export path, so it stays editable.
  const showKindEditor = showAgentConfig || showWorkflow;
  const headerHiddenFields = showAgentConfig
    ? ['resource', 'type', DATA_FORMAT_FIELD]
    : (showWorkflow ? ['type', DATA_FORMAT_FIELD] : []);

  const isChanged = JSON.stringify(doc.data) !== JSON.stringify(shownDoc);
  return (
    <>
      <Stack direction={'row'} alignItems={'center'} justifyItems={'center'}>
        <Typography variant='h6' sx={{ marginRight: 1 }}>
          {doc.isConfig ? doc.project.name + ' config' : doc.name}
        </Typography>
        <DocViewSelector
          docView={docView}
          setDocView={setDocView}
          enabled={{ [TabKind.Agent]: isAgent, [TabKind.Workflow]: isWorkflow }}
        />
        <DetailsVisibilityToggle
          value={detailsVisibility}
          onChange={setDetailsVisibility}
        />
        {isWorkflow && (
          <RunWorkflowButton
            projectName={doc.project.name}
            workflowName={(shownDoc.desc as WorkflowDesc).workflowName ?? doc.name}
            disabled={isChanged}
            disabledReason="Save changes before running"
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
      </Stack>
      {showHeader && (
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
      )}
      {showTree && (
        <SimpleTreeView
          defaultExpandedItems={[keyForDetailsViewItem(DESC_FIELD), keyForDetailsViewItem('resource')]}
          // Rows reserve extra space below for the "required" helper text, which
          // centers the chevron a bit low; nudge it up to line up with the name.
          sx={{ '& .MuiTreeItem-iconContainer': { transform: 'translateY(-3px)' } }}
        >
          {reorderEntries(Object.entries(shownDoc), [DESC_FIELD, 'resource']).map(([k, v]) => {
          if (FORBIDDEN_FIELDS.includes(k)) {
            return null;
          }
          if (headerHiddenFields.includes(k)) {
            return null;
          }
          const hideOnDesc = showFormulated && k === DESC_FIELD;
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
      )}
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
              disabled={isChanged}
              disabledReason="Save changes before running"
              sx={{ bgcolor: 'background.paper', boxShadow: 1, p: 0.5 }}
            />
          }
        />
      )}
    </>
  );
}
