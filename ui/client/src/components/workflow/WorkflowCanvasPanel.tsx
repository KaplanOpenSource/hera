import { Box, Typography } from '@mui/material';
import { DocumentObj, ProjectObj } from '../../objects/ProjectObj';
import { WorkflowDesc } from '../../shared/types';
import { useShownDoc } from '../details/useShownDoc';
import { RunWorkflowButton } from './RunWorkflowButton';
import { WorkflowEditor } from './WorkflowEditor';

// One workflow document's canvas, as its own dock tab. It edits the same shown
// document as the details tab, so the two stay in step.
export const WorkflowCanvasPanel = ({
  project,
  docid,
}: {
  project: ProjectObj,
  docid: string,
}) => {
  const doc = project.allDocuments.find(d => d.docid === docid);

  return doc
    ? <WorkflowCanvas project={project} doc={doc} />
    : <Typography sx={{ p: 2 }} color="text.secondary">Document not found.</Typography>;
};

const WorkflowCanvas = ({
  project,
  doc,
}: {
  project: ProjectObj,
  doc: DocumentObj,
}) => {
  const { shownDoc, setShownDoc, saveShownDoc } = useShownDoc(doc);
  const isChanged = JSON.stringify(doc.data) !== JSON.stringify(shownDoc);

  return (
    <Box sx={{ height: '100%', display: 'flex', flexDirection: 'column', minHeight: 0 }}>
      <WorkflowEditor
        workflow={(shownDoc.desc as WorkflowDesc).workflow}
        setWorkflow={newVal => setShownDoc({ ...shownDoc, desc: { ...shownDoc.desc, workflow: newVal } as WorkflowDesc })}
        actionButtons={
          <RunWorkflowButton
            projectName={project.name}
            workflowName={(shownDoc.desc as WorkflowDesc).workflowName ?? doc.name}
            doc={shownDoc}
            isChanged={isChanged}
            save={saveShownDoc}
            sx={{ bgcolor: 'background.paper', boxShadow: 1, p: 0.25, '& .MuiSvgIcon-root': { fontSize: 16 } }}
          />
        }
      />
    </Box>
  );
};
