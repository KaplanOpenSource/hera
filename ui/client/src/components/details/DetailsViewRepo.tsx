import { Edit } from "@mui/icons-material";
import { Stack } from "@mui/material";
import { SimpleTreeView, TreeItem } from "@mui/x-tree-view";
import { useEffect, useState } from "react";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { JsonTreeNode } from "../../elements/JsonTreeView";
import { execPython } from "../../io/execPython";
import { idRepoId } from "../../shared/idDocId";
import { RepoTreeAddButton } from "./RepoTreeAddButton";
import { TextFieldWithApply } from "./TextFieldWithApply";

export const DetailsViewRepo = ({
  repoPath,
}: {
  repoPath: string,
}) => {
  const [tree, setTree] = useState<any>(undefined);
  const [repoStr, setRepoStr] = useState<string>('');
  const [showStr, setShowStr] = useState<boolean>(false);

  useEffect(() => {
    (async () => {
      const { data } = await execPython(`
import json
with open('${repoPath}', 'r') as fjson:
  data = json.load(fjson)
result = {"json": data}
            `);
      console.log(data)
      if (data?.json) {
        setTree(data.json);
      } else {
        setTree(undefined)
      }
    })()
  }, [repoPath])

  useEffect(() => {
    const newText = JSON.stringify(tree, undefined, 2);
    if (newText !== repoStr) {
      setRepoStr(JSON.stringify(tree, undefined, 2))
    }
  }, [tree])

  return (
    <>
      <SimpleTreeView
        defaultExpandedItems={[idRepoId(repoPath), 'root.root']}
      >
        <TreeItem
          key={idRepoId(repoPath)}
          itemId={idRepoId(repoPath)}
          label={(
            <Stack direction={'row'} justifyItems={'center'} alignItems={'center'}>
              {repoPath}
              <RepoTreeAddButton
                tree={tree}
              />
              <ButtonTooltip
                title={'Edit Repository as string Json'}
                onClick={() => setShowStr(!showStr)}
              >
                <Edit />
              </ButtonTooltip>
            </Stack>
          )}
        >
          {tree === undefined ? null : (
            <JsonTreeNode
              label="root"
              value={tree}
              setData={val => {
                setTree(val)
              }}
              parentKey="root"
            />
          )}
        </TreeItem>
      </SimpleTreeView>
      {showStr
        ? (
          <TextFieldWithApply
            value={repoStr}
            setValue={v => setRepoStr(v)}
            onApply={() => {
              try {
                const t = JSON.parse(repoStr);
                setTree(t);
              } catch {
              }
            }}
            textTitle="Repository Json (as string)"
            buttonTitle={'Apply string Json to tree'}
          />
        )
        : null}
    </>
  )
}