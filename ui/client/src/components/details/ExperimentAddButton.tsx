import { RepoTreeAddButton } from "./RepoTreeAddButton";

export const ExperimentAddButton = ({
  tree,
  treePath,
}: {
  tree: any,
  treePath: string[],
}) => {
  const [section, expName] = treePath;
  const expValue = tree?.[section]?.[expName];
  return (
    <RepoTreeAddButton
      tree={{ [section]: { [expName]: expValue } }}
      title={`Add experiment '${expName}' to project`}
    />
  );
};
