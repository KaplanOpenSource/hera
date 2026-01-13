import { AccountTree, Delete } from "@mui/icons-material";
import { Stack } from "@mui/material";
import { TreeItem } from "@mui/x-tree-view";
import { useState } from "react";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { JsonTreeNode } from "../../elements/JsonTreeView";
import { execPython } from "../../io/execPython";

export const RepoTreeOne = ({
  repoPath,
  setRepoPath,
}: {
  repoPath: string,
  setRepoPath: (v: string | undefined) => void,
}) => {
  const [tree, setTree] = useState<any>(undefined);

  const showRepo = async () => {
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
  }

  return (
    <TreeItem
      key={'__repo_*_' + repoPath}
      itemId={'__repo_*_' + repoPath}
      label={(
        <Stack direction={'row'} justifyItems={'center'} alignItems={'center'}>
          {repoPath}
          <ButtonTooltip
            title={'Remove repository'}
            onClick={() => setRepoPath(undefined)}
          >
            <Delete />
          </ButtonTooltip>
          <ButtonTooltip
            title={'Show repository'}
            onClick={() => showRepo()}
          >
            <AccountTree />
          </ButtonTooltip>
          {/* {JSON.stringify(tree)} */}
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
  )
}