import { TreeItem } from "@mui/x-tree-view"
import { ProjectDocument, ProjectEntire } from "@shared/types"

export const ProjectDocumentItem = ({
  project,
  document,
}: {
  project: ProjectEntire,
  document: ProjectDocument,
}) => {
  const id = `document${document.desc.docid}`;
  return (
    <TreeItem
      key={id} itemId={id} label={`Document: ${document.desc.datasourceName}`}
    >
      <TreeItem
        key={`${id}-version`} itemId={`${id}-version`} label={`Version: ${document.desc.version.join('.')}`}
      >
      </TreeItem>

    </TreeItem>
  )
}