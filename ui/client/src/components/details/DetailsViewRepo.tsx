import { Box, Stack, TextField } from "@mui/material";
import { SimpleTreeView, TreeItem } from "@mui/x-tree-view";
import { useEffect, useState } from "react";
import { JsonTreeNode } from "../../elements/JsonTreeView";
import { execPython } from "../../io/execPython";
import { idRepoId } from "../../shared/idDocId";
import { RepoTreeAddButton } from "./RepoTreeAddButton";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { Check, Edit } from "@mui/icons-material";

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
          <Box position={'relative'}>
            <TextField
              fullWidth
              label="Repository Json (as string)"
              value={repoStr}
              onClick={e => e.stopPropagation()}
              onChange={(e) => {
                setRepoStr(e.target.value);
              }}
              size="small"
              rows={15}
              multiline={true}
            />
            <ButtonTooltip
              sx={{
                position: 'absolute',
                top: '8px',
                right: '20px',
                zIndex: 1
              }}
              title={'Apply string Json to tree'}
              onClick={() => {
                try {
                  const t = JSON.parse(repoStr);
                  setTree(t);
                } catch {
                }
              }}
            >
              <Check />
            </ButtonTooltip>
          </Box>
        )
        : null}
    </>
  )
}