import { Delete } from "@mui/icons-material";
import { Box, Stack, Typography } from "@mui/material";
import { TreeItem } from "@mui/x-tree-view";
import { ProjectDocument, ProjectEntire } from "@shared/types";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { useConfirm } from "../../elements/useConfirm";
import { fetchProjectDetails } from "../../io/FetchProjects";
import { execPython } from "../../io/execPython";
import { idDocId } from "../../shared/idDocId";

export const ProjectDocumentItem = ({
  project,
  document,
}: {
  project: ProjectEntire,
  document: ProjectDocument,
}) => {
  const { confirmOpen, ConfirmDialog } = useConfirm()

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
        <Stack direction='column' justifyContent="start" alignItems=''>
          <Stack direction='row' spacing={1} justifyContent="start" alignItems='center'>
            <Typography>
              {name}
            </Typography>
            <ButtonTooltip
              title={isProjectConfig ? 'Project Config is deleted only with project' : 'Delete Document'}
              disabled={isProjectConfig}
              onClick={async () => {
                if ((await confirmOpen({ title: `Are you sure you want to delete ${name}?` })).confirmed) {
                  deleteDocument()
                }
              }}
            >
              <Delete />
              {ConfirmDialog}
            </ButtonTooltip>
          </Stack>
          {document?.resource && typeof (document?.resource) == 'string' && (
            <Typography sx={{fontSize: 10}}>
              resource: {document.resource.substring(0, 80)}
            </Typography>
          )}
        </Stack>
      }
    >
    </TreeItem>
  )
}
