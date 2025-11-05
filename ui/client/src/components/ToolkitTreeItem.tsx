import { TreeItem } from "@mui/x-tree-view";
import { ProjectDocumentItem } from "./ProjectDocumentItem";
import { ProjectEntire, Toolkit } from "@shared/types";

export const ToolkitTreeItem = ({
  project,
  toolkit,
}: {
  toolkit: Toolkit | null,
  project: ProjectEntire,
}) => {
  const documentsForToolkit = (toolkitName: string) => {
    return project?.documents.filter(d => d?.desc?.toolkit === toolkitName) || [];
  }

  const toolkitName = toolkit?.toolkit || '';
  const docs = documentsForToolkit(toolkitName);
  return !docs.length ? null : (
    <TreeItem key={toolkitName} itemId={toolkitName} label={toolkitName}>
      {/* <Typography>
          {cls}
        </Typography>
        {desc ?? (
          <Typography>
            {desc}
          </Typography>
        )} */}
      {docs.map(d => {
        return (
          <ProjectDocumentItem
            key={`proj${project.name}_doc${d._id.$oid}`} project={project} document={d}
          >
          </ProjectDocumentItem>
        )
      })}
    </TreeItem>
  )
}