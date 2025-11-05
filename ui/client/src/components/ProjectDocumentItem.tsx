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
  const name = document?.desc?.datasourceName || document?.type || document._cls;
  // console.log(document.desc.version)
  return (
    <TreeItem
      key={id} itemId={id}
      label={name}
    >
      <TreeItem
        key={id + '-details'} itemId={id + '-details'} label={
          <Box>
            {document.desc?.version && (
              <Typography>
                Version: {(document.desc?.version || []).join('.')}
              </Typography>
            )}
            <Typography>
              Type: {document.type}
            </Typography>
            {document?.resource && typeof (document?.resource) == 'string' && (
              <Typography>
                resource: {document.resource}
              </Typography>
            )}
            {document.desc.toolkit && (
              <Typography>
                toolkit: {document.desc.toolkit}
              </Typography>
            )}
          </Box>
        }
      >
      </TreeItem>
    </TreeItem>
  )
}