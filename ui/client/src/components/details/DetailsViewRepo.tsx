import { Delete, Edit } from "@mui/icons-material";
import { Stack, Tooltip, Typography } from "@mui/material";
import { SimpleTreeView, TreeItem } from "@mui/x-tree-view";
import { useEffect, useState } from "react";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { fetchPython } from "../../io/fetchPython";
import { idRepoId, TEMP_REPO_NAME } from "../../shared/idDocId";
import { RepoTreeAddButton } from "./RepoTreeAddButton";
import { TextFieldWithApply } from "./TextFieldWithApply";
import { JsonTreeNode } from "../../elements/JsonTreeNode";

export const DetailsViewRepo = ({
  repoPath,
}: {
  repoPath: string,
}) => {
  const isTempRepo = repoPath.includes(TEMP_REPO_NAME);
  const [tree, setTree] = useState<any>(undefined);
  const [repoStr, setRepoStr] = useState<string>('');
  const [showStr, setShowStr] = useState<boolean>(isTempRepo);

  useEffect(() => {
    (async () => {
      if (!isTempRepo) {
        const { data } = await fetchPython({
          results: ['jsonData'],
          label: 'repo JSON',
          code: `
import json
with open('${repoPath}', 'r') as fjson:
  jsonData = json.load(fjson)
`,
        });
        if (data?.jsonData) {
          setTree(data.jsonData);
        } else {
          setTree(undefined)
        }
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
              <Tooltip
                title={<>
                  Json is in Folder:<br />
                  {repoPath.split('/').slice(0, -1).join('/')}
                </>}
              >
                <Typography marginRight={1}>
                  {repoPath.split('/').pop()}
                </Typography>
              </Tooltip>
              <RepoTreeAddButton
                tree={tree}
              />
              <ButtonTooltip
                title={'Edit Repository as Json text'}
                onClick={() => setShowStr(!showStr)}
              >
                <Edit />
              </ButtonTooltip>
            </Stack>
          )}
        >
          {typeof tree === "object" && tree !== null
            ? (
              <>
                {Object.entries(tree).map(([key, val]) => (
                  <JsonTreeNode
                    key={key}
                    label={key}
                    parentKey={`root.${key}`}
                    value={val as any}
                    rootLabelComponents={(
                      <ButtonTooltip
                        onClick={() => setTree(Object.fromEntries(Object.entries(tree).filter(([k, _]) => k !== key)))}
                        title={'Remove Data Source from repository (locally, not affecting the file)'}
                      >
                        <Delete fontSize="inherit" />
                      </ButtonTooltip>
                    )}
                  />
                ))}
              </>
            )
            : null}
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
            textTitle="Repository Json (as text)"
            buttonTitle={'Apply Json text to tree (locally, not affecting the file)'}
          />
        )
        : null}
    </>
  )
}