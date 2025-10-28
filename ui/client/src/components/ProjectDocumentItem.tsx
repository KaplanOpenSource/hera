import { Box, Typography } from "@mui/material";
import { TreeItem } from "@mui/x-tree-view"
import { ProjectDocument, ProjectEntire } from "@shared/types"

export const ProjectDocumentItem = ({
  project,
  document,
}: {
  project: ProjectEntire,
  document: ProjectDocument,
}) => {
  const id = `document${document?.desc?.docid}`;
  return (
    <TreeItem
      key={id} itemId={id} label={`Document: ${document.desc.datasourceName}`}
    >
      <TreeItem
        key={id + '-details'} itemId={id + '-details'} label={
          <>
            <Typography>
              Version: {(document.desc?.version || []).join('.')}
            </Typography>
            <Typography>
              Type: {document.type}
            </Typography>
            <Typography>
              resource: {document.resource}
            </Typography>
            <Typography>
              toolkit: {document.desc.toolkit}
            </Typography>
          </>
        }
      >
      </TreeItem>
    </TreeItem>
  )
}