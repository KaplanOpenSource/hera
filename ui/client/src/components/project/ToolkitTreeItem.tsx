import { TreeItem } from "@mui/x-tree-view";
import { ProjectDocumentItem } from "./ProjectDocumentItem";
import { ProjectEntire, Toolkit } from "@shared/types";
import { useProjectStore } from "../../stores/useProjectStore";
import { Add } from "@mui/icons-material";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { Stack, Typography } from "@mui/material";
import { AddDocumentButton } from "./AddDocumentButton";

export const ToolkitTreeItem = ({
  project,
  toolkit,
  showEmpty = false,
}: {
  toolkit: Toolkit | undefined,
  project: ProjectEntire,
  showEmpty?: boolean,
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
      {docs.map(d => {
        return (
          <ProjectDocumentItem
            key={`proj${project.name}_doc${d._id.$oid}`}
            project={project}
            document={d}
          />
        )
      })}
    </TreeItem>
  )
}