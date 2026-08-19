import { Stack, Typography } from "@mui/material";
import { TreeItem } from "@mui/x-tree-view";
import { ProjectDocument, ProjectEntire } from "@shared/types";
import { idDocId } from "../../shared/idDocId";
import { TabKind } from "../../shared/tabKind";
import { useViewSettingsStore } from "../../stores/useViewSettingsStore";
import { DocumentKindIcon } from "./DocumentKindIcon";

// The per-document delete button now lives in the document interface
// (see DeleteDocumentButton in the details view), not in this tree item.
export const ProjectDocumentItem = ({
  document,
  kind,
}: {
  project: ProjectEntire,
  document: ProjectDocument,
  kind: TabKind,
  onDocumentDeleted?: () => void,
}) => {
  const { viewSettings } = useViewSettingsStore();

  const id = idDocId(document?._id.$oid);
  const name = document?.desc?.datasourceName || document?.type || document._cls;

  return (
    <TreeItem
      key={id} itemId={id}
      label={
        <Stack direction='column' justifyContent="start" alignItems=''>
          <Stack direction='row' spacing={1} justifyContent="start" alignItems='center'>
            <DocumentKindIcon kind={kind} fontSize="small" sx={{ color: 'text.secondary' }} />
            <Typography>
              {name}
            </Typography>
          </Stack>
          {viewSettings.showDocumentPreview && document?.resource && typeof (document?.resource) == 'string' && (
            <Typography sx={{ fontSize: 10 }}>
              resource: {document.resource.substring(0, 80)}
            </Typography>
          )}
        </Stack>
      }
    >
    </TreeItem>
  )
}
