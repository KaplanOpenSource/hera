import { SECTION_EXPERIMENT } from "./RepoJsonMerger";
import { RepoTreeAddButton } from "./RepoTreeAddButton";

export const ToolkitAddButton = ({
  tree,
  toolkitName,
}: {
  tree: any,
  toolkitName: string,
}) => {
  const sections = tree?.[toolkitName];
  if (!sections || typeof sections !== 'object') return null;
  const filtered = Object.fromEntries(
    Object.entries(sections).filter(([key]) => key.toLowerCase() !== SECTION_EXPERIMENT)
  );
  return (
    <RepoTreeAddButton
      tree={{ [toolkitName]: filtered }}
      title={`Add datasources of toolkit '${toolkitName}' to project`}
    />
  );
};
