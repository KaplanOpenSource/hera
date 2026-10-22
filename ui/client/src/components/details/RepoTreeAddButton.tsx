import { DriveFolderUpload } from "@mui/icons-material";
import { BooleanProperty } from "../../elements/BooleanProperty";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { useDialog } from "../../elements/useDialog";
import { loadRepositoryJsonIntoProject } from "../../io/loadRepositoryIntoProject";
import { useProjectStore } from "../../stores/useProjectStore";

const filterTree = (obj: any, hidden: Set<string>, prefix: string[] = []): any => {
  if (typeof obj !== 'object' || obj === null || Array.isArray(obj)) return obj;
  const result: { [key: string]: any } = {};
  for (const [key, val] of Object.entries(obj)) {
    const pathKey = [...prefix, key].join('/');
    if (!hidden.has(pathKey)) {
      result[key] = filterTree(val, hidden, [...prefix, key]);
    }
  }
  return result;
};

export const RepoTreeAddButton = ({
  tree,
  hiddenPaths,
  title = 'Add repository of data sources to project',
}: {
  tree: any,
  hiddenPaths?: Set<string>,
  title?: string,
}) => {
  type AddRepoArgs = {
    overwrite: boolean;
  };

  const { openDialog, DialogComponent } = useDialog<AddRepoArgs>();
  const { currProject } = useProjectStore();

  const baseDir = currProject?.documents.find(x => x.type === currProject.name + '__config__')?.desc.filesDirectory || '';

  const addRepo = async (params: AddRepoArgs) => {
    // filterTree is about what the user hid in the browser, so it stays here.
    const visibleTree = hiddenPaths?.size ? filterTree(tree, hiddenPaths) : tree;
    await loadRepositoryJsonIntoProject(visibleTree, baseDir, params.overwrite);
  }

  return (
    <ButtonTooltip
      color="primary"
      title={title}
      disabled={tree === undefined}
      onClick={async () => {
        const result = await openDialog({
          title,
          initialValues: { overwrite: true },
          render:
            ({ values, setValues }) => (
              <BooleanProperty
                label="Overwrite"
                value={values.overwrite}
                setValue={(v) => setValues({ ...values, overwrite: v })}
              />
            )
        });

        if (result.confirmed && result.values) {
          await addRepo(result.values)
        }
      }}
    >
      <DriveFolderUpload />
      {DialogComponent}
    </ButtonTooltip>
  )
}