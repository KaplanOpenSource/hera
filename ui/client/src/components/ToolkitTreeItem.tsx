import { TreeItem } from "@mui/x-tree-view";
import { ProjectDocumentItem } from "./ProjectDocumentItem";
import { ProjectEntire, Toolkit } from "@shared/types";
import { useProjectStore } from "../stores/useProjectStore";

export const ToolkitTreeItem = ({
  project,
  toolkit,
}: {
  toolkit: Toolkit | undefined,
  project: ProjectEntire,
}) => {
  const { toolkits } = useProjectStore();

  const documentsForToolkit = (toolkitName: string) => {
    return project?.documents.filter(d => d?.desc?.toolkit === toolkitName) || [];
  }

  const documentsWithoutToolkit = () => {
    return project?.documents.filter(d => {
      const found = toolkits.find(({ toolkit }) => toolkit === d?.desc?.toolkit);
      return found === undefined;
    });
  }

  const toolkitName = toolkit ? toolkit.toolkit : 'no-toolkit';
  const toolkitLabel = toolkit ? toolkit.toolkit : 'No Toolkit Documents';
  const docs = toolkit ? documentsForToolkit(toolkitName) : documentsWithoutToolkit();
  return !docs.length ? null : (
    <TreeItem key={toolkitName} itemId={toolkitName} label={toolkitLabel}>
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