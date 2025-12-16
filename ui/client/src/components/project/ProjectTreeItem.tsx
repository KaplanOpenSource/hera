import { TreeItem } from '@mui/x-tree-view/TreeItem';
import type { Project } from '../../shared/types';

export const ProjectTreeItem = ({ project }: { project: Project }) => {
  const hasDetails = project.documents && project.toolkitCount !== undefined;

  return (
    <TreeItem key={project.id} itemId={project.id} label={project.name}>
      {hasDetails
        ? [
          <TreeItem key={`${project.id}-docs`} itemId={`${project.id}-docs`} label="Documents">
            <TreeItem
              key={`${project.id}-docs-sim`}
              itemId={`${project.id}-docs-sim`}
              label={`Sim: ${project.documents!.sim}`}
            />
            <TreeItem
              key={`${project.id}-docs-measure`}
              itemId={`${project.id}-docs-measure`}
              label={`Measure: ${project.documents!.measure}`}
            />
            <TreeItem
              key={`${project.id}-docs-cache`}
              itemId={`${project.id}-docs-cache`}
              label={`Cache: ${project.documents!.cache}`}
            />
          </TreeItem>,
          <TreeItem
            key={`${project.id}-toolkits`}
            itemId={`${project.id}-toolkits`}
            label={`Toolkits: ${project.toolkitCount}`}
          />,
        ]
        : <TreeItem key={`${project.id}-loading`} itemId={`${project.id}-loading`} label="Loading..." disabled />}
    </TreeItem>
  );
};
