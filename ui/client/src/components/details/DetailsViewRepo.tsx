import { Stack } from "@mui/material";
import { SimpleTreeView, TreeItem } from "@mui/x-tree-view";
import { useEffect, useState } from "react";
import { JsonTreeNode } from "../../elements/JsonTreeView";
import { execPython } from "../../io/execPython";
import { idRepoId } from "../../shared/idDocId";
import { RepoTreeAddButton } from "./RepoTreeAddButton";

export const DetailsViewRepo = ({
  repoPath,
}: {
  repoPath: string,
}) => {
  const [tree, setTree] = useState<any>(undefined);

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

  return (
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
  )
}