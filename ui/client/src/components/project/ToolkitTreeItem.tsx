import { Stack, Typography } from "@mui/material";
import { TreeItem } from "@mui/x-tree-view";
import { Toolkit } from "@shared/types";
import { ProjectObj } from "../../objects/ProjectObj";
import { useProjectStore } from "../../stores/useProjectStore";
import { AddDocumentButton } from "./AddDocumentButton";
import { ProjectDocumentItem } from "./ProjectDocumentItem";

export const ToolkitTreeItem = ({
  project,
  toolkit,
  showEmpty = false,
  showDocumentPreview = true,
}: {
  toolkit: Toolkit | undefined,
  project: ProjectObj,
  showEmpty?: boolean,
  showDocumentPreview?: boolean,
}) => {
  const { toolkits } = useProjectStore();

  const documentsForToolkit = (toolkitName: string) => {
    return project?.documents.filter(d => d.toolkit === toolkitName) || [];
  }

  const documentsWithoutToolkit = () => {
    return project?.documents.filter(d => {
      const found = toolkits.find(({ toolkit }) => toolkit === d.toolkit);
      return found === undefined;
    });
  }

  const toolkitName = toolkit ? toolkit.toolkit : 'no-toolkit';
  const toolkitLabel = toolkit ? toolkit.toolkit : 'No Toolkit Documents';
  const docs = toolkit ? documentsForToolkit(toolkitName) : documentsWithoutToolkit();
  return (!docs.length && !showEmpty) ? null : (
    <TreeItem key={toolkitName} itemId={toolkitName}
      label={
        <Stack direction='row' spacing={1} justifyContent="start" alignItems='center'>
          <Typography>
            {toolkitLabel}
          </Typography>
          <AddDocumentButton
            toolkit={toolkit}
          />
        </Stack>
      }
    >
      {/* <Typography>
          {cls}
        </Typography>
        {desc ?? (
          <Typography>
            {desc}
          </Typography>
        )} */}
      {docs.map(document => {
        return (
          <ProjectDocumentItem
            key={`proj${project.name}_doc${document.docid}`}
            project={project.data}
            document={document.data}
            showDocumentPreview={showDocumentPreview}
          />
        )
      })}
    </TreeItem>
  )
}