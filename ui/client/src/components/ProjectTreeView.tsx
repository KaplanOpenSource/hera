import { Paper, Typography } from '@mui/material';
import { SimpleTreeView } from '@mui/x-tree-view/SimpleTreeView';
import type { Project, ProjectEntire } from '../shared/types';
import { ProjectTreeItem } from './ProjectTreeItem';
import { TreeItem } from '@mui/x-tree-view';
import { ProjectDocumentItem } from './ProjectDocumentItem';
import { useState } from 'react';
import { FetchPython } from '../io/FetchPython';

export const ProjectTreeView = ({
  project,
  // onProjectSelect,
  // onProjectExpand,
}: {
  project: ProjectEntire | null;
  // onProjectSelect: (project: Project) => void;
  // onProjectExpand: (projectId: string) => void;
}) => {
  // const handleItemSelection = (event: React.SyntheticEvent | null, itemId: string | null) => {
  //   if (itemId) {
  //     const selectedProject = projects.find((p) => p.id === itemId);
  //     if (selectedProject) {
  //       onProjectSelect(selectedProject);
  //     }
  //   }
  // };

  // const handleItemExpansion = (event: React.SyntheticEvent | null, itemIds: string[]) => {
  //   const lastExpanded = itemIds[0];
  //   if (lastExpanded) {
  //     onProjectExpand(lastExpanded);
  //   }
  // };

  const [toolkits, setToolKits] = useState<[string, any][]>([]);

  console.log(toolkits)
  return (
    <Paper sx={{ p: 2 }}>
      <FetchPython
        code={`
from hera import toolkitHome
import json
result = json.dumps(toolkitHome._toolkits,indent=4)
        `}
        onSuccess={(data) => {
          setToolKits(Object.entries(JSON.parse(data)))
        }}
      />
      <SimpleTreeView
      >
        {!project
          ? (
            <TreeItem key={`no-project`} itemId={`no-project`} label="No project loaded" disabled />
          ) : (
            <TreeItem key={`project-documents`} itemId={`project-documents`} label={`Project ${project.name}`}>
              {project.documents.map(d => (
                <ProjectDocumentItem
                  key={`proj${project.name}_doc${d.desc.docid}`} project={project} document={d}
                >
                </ProjectDocumentItem>
              ))}
            </TreeItem>
          )
        }
        <TreeItem key={`toolkits`} itemId={`toolkits`} label={`Toolkits`}>
          {toolkits.map(([name, { cls, desc }]) => (
            <TreeItem key={name} itemId={name} label={name}>
              <Typography>
                {cls}
              </Typography>
              {desc ?? (
                <Typography>
                  {desc}
                </Typography>
              )}
            </TreeItem>
          ))}
        </TreeItem>
      </SimpleTreeView>
    </Paper>
  );
};
