import { Stack, Typography } from "@mui/material";
import { TreeItem } from "@mui/x-tree-view";
import { Toolkit } from "@shared/types";
import { ProjectObj } from "../../objects/ProjectObj";
import { useProjectStore } from "../../stores/useProjectStore";
import { AddDocumentButton } from "./AddDocumentButton";
import { ProjectDocumentItem } from "./ProjectDocumentItem";
import { compareJsons, filterAndSortByGroups, sortByDistinctValues } from "../../utils/compareJsons";

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

  const toolkitName = toolkit ? toolkit.toolkit : 'no-toolkit';
  const toolkitLabel = toolkit ? toolkit.toolkit : 'No Toolkit Documents';
  const docs = project?.documents.filter(d => {
    if (toolkit) {
      return d.toolkit === toolkitName;
    }
    return !toolkits.some(t => t.toolkit === d.toolkit);
  }) || [];

  console.log('compareJsons:\n', filterAndSortByGroups(compareJsons(docs.map(x => x.data.desc)), 2));
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