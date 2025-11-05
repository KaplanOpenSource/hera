import { Box, Stack, Typography } from "@mui/material";
import { TreeItem } from "@mui/x-tree-view"
import { ProjectDocument, ProjectEntire } from "@shared/types"
import { ButtonTooltip } from "../elements/ButtonTooltip";
import { Delete } from "@mui/icons-material";
import { fetchProjectDetails } from "../io/FetchProjects";
import { execPython } from "../io/execPython";

export const ProjectDocumentItem = ({
  project,
  document,
}: {
  project: ProjectEntire,
  document: ProjectDocument,
}) => {
  const id = `document${document?._id.$oid}`;
  const name = document?.desc?.datasourceName || document?.type || document._cls;

  const deleteDocument = async () => {
    const { problem } = await execPython(`
from hera.datalayer import All
All.deleteDocumentByID('${document?._id.$oid}')
        `)
    if (problem) {
      return;
    }
    await fetchProjectDetails(project.name);
  }

  return (
    <TreeItem
      key={id} itemId={id}
      label={
        <Stack direction='row' spacing={1} justifyContent="start" alignItems='center'>
          <Typography>
            {name}
          </Typography>
          <ButtonTooltip
            title={'Delete Document'}
            onClick={() => {
              if (confirm(`Are you sure you want to delete ${name}?`)) {
                deleteDocument()
              }
            }}
          >
            <Delete />
          </ButtonTooltip>
        </Stack>

      }
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