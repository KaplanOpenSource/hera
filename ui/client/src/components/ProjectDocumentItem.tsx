import { Box, Stack, Typography } from "@mui/material";
import { TreeItem } from "@mui/x-tree-view"
import { ProjectDocument, ProjectEntire } from "@shared/types"
import { ButtonTooltip } from "../elements/ButtonTooltip";
import { Delete } from "@mui/icons-material";
import { fetchProjectDetails } from "../io/FetchProjects";
import { execPython } from "../io/execPython";
import { idDocId } from "../shared/idDocId";

export const ProjectDocumentItem = ({
  project,
  document,
}: {
  project: ProjectEntire,
  document: ProjectDocument,
}) => {
  const id = idDocId(document?._id.$oid);
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

  const isProjectConfig = document.type === project.name + '__config__'

  return (
    <TreeItem
      key={id} itemId={id}
      label={
        <Stack direction='row' spacing={1} justifyContent="start" alignItems='center'>
          <Typography>
            {name}
          </Typography>
          <ButtonTooltip
            title={isProjectConfig ? 'Project Config is deleted only with project' : 'Delete Document'}
            disabled={isProjectConfig}
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