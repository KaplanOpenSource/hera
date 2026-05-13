import { RepoTreeAddButton } from "./RepoTreeAddButton";

export const DatasourceAddButton = ({
  tree,
  treePath,
}: {
  tree: any,
  treePath: string[],
}) => {
  const [toolkit, dsSection, dsName] = treePath;
  const dsValue = tree?.[toolkit]?.[dsSection]?.[dsName];
  const singleDsTree = { [toolkit]: { [dsSection]: { [dsName]: dsValue } } };
  return (
    <RepoTreeAddButton
      tree={singleDsTree}
      title={`Add datasource '${dsName}' to project`}
    />
  );
};
