import { Stack, Typography } from "@mui/material";
import { TreeItem } from "@mui/x-tree-view";
import { Toolkit } from "@shared/types";
import { ProjectObj } from "../../objects/ProjectObj";
import { useProjectStore } from "../../stores/useProjectStore";
import { AddDocumentButton } from "./AddDocumentButton";
import { DocumentSplitGroup } from "./DocumentSplitGroup";
import { ViewSettingsType } from "./ProjectTreeView";

export const ToolkitTreeItem = ({
  project,
  toolkit,
  viewSettings,
}: {
  toolkit: Toolkit | undefined;
  project: ProjectObj;
  viewSettings: ViewSettingsType;
}) => {
  const { toolkits } = useProjectStore();

  const toolkitName = toolkit ? toolkit.toolkit : "no-toolkit";
  const toolkitLabel = toolkit ? toolkit.toolkit : "No Toolkit Documents";
  const docs = (project?.documents || []).filter((d) => {
    if (toolkit) {
      return d.toolkit === toolkitName;
    }
    return !toolkits.some((t) => t.toolkit === d.toolkit);
  });

  if (docs.length === 0 && !viewSettings.showEmptyToolkits) {
    return null;
  }

  return (
    <TreeItem
      key={toolkitName}
      itemId={toolkitName}
      label={
        <Stack direction="row" spacing={1} justifyContent="start" alignItems="center">
          <Typography>{toolkitLabel}</Typography>
          <AddDocumentButton toolkit={toolkit} />
        </Stack>
      }
    >
      <DocumentSplitGroup
        docs={docs}
        project={project}
        showDocumentPreview={viewSettings.showDocumentPreview}
        depth={viewSettings.maxDepth}
        minGroupSize={viewSettings.minGroupSize}
      />
    </TreeItem>
  );
};
