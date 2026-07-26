import { ContentCopy } from "@mui/icons-material";
import { Stack, TextField } from "@mui/material";
import { useNavigate } from "react-router-dom";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { useDialog } from "../../elements/useDialog";
import { fetchPython } from "../../io/fetchPython";
import { ProjectName } from "../../shared/types";
import { useProjectStore } from "../../stores/useProjectStore";

interface CopyProjectValues {
  newName: string;
  newDir: string;
}

export const CopyProjectButton = ({
  onCopied,
}: {
  onCopied?: () => void,
}) => {
  const { openDialog, DialogComponent } = useDialog<CopyProjectValues>();
  const { currProjectName, setProjectNames } = useProjectStore();
  const navigate = useNavigate();

  // Default the copy into a sibling folder of the source, named after the copy.
  const srcFilesDir = useProjectStore.getState().getProject()?.configDocument?.data.desc.filesDirectory ?? '';
  const defaultName = `${currProjectName}_copy`;
  const parentDir = srcFilesDir.replace(/\/+$/, '').split('/').slice(0, -1).join('/');
  const defaultDir = parentDir
    ? `${parentDir}/${defaultName}`
    : (srcFilesDir ? `${srcFilesDir}/${defaultName}` : '');

  const copyProject = async (newName: string, newDir: string) => {
    const dirExpr = newDir !== ''
      ? `'${newDir}'`
      : `os.path.join(os.getcwd(), 'projects', '${newName}')`;
    const { data } = await fetchPython({
      results: ['projectNames'],
      label: `copy project ${currProjectName} to ${newName}`,
      code: `
import os, shutil
from hera.datalayer import All
from hera.datalayer.project import getProjectList, createProjectDirectory, Project

srcName = '${currProjectName}'
newName = '${newName}'
if newName in getProjectList():
    raise ValueError(f"A project named {newName} already exists")

srcDir = os.path.abspath(Project(projectName=srcName).filesDirectory or '')
newDir = os.path.abspath(${dirExpr})

# Create the destination project (writes caseConfiguration.json and the config document).
createProjectDirectory(outputPath=newDir, projectName=newName)
newProj = Project(projectName=newName, filesDirectory=newDir)

configType = srcName + '__config__'
for doc in All.getDocuments(projectName=srcName):
    d = doc.asDict(with_id=False)
    if d.get('type') == configType:
        continue  # the destination already has its own config document
    res = d.get('resource')
    # Copy single on-disk files that live inside the source files directory, preserving their relative path.
    # Directory resources are left pointing at the original folder for now (not copied).
    if isinstance(res, str) and res and srcDir:
        resAbs = os.path.abspath(res)
        if (resAbs == srcDir or resAbs.startswith(srcDir + os.sep)) and os.path.isfile(resAbs):
            newRes = os.path.join(newDir, os.path.relpath(resAbs, srcDir))
            os.makedirs(os.path.dirname(newRes), exist_ok=True)
            shutil.copy2(resAbs, newRes)
            d['resource'] = newRes
    newProj.addDocumentFromDict(d)

projectNames = [{"name": proj} for proj in getProjectList()]
`,
    });
    if (!data) {
      return;
    }
    setProjectNames((data.projectNames || []) as ProjectName[]);
    navigate('/' + encodeURIComponent(newName));
    onCopied?.();
  }

  return (<>
    <ButtonTooltip
      button
      title='Copy project'
      onClick={async () => {
        const { confirmed, values } = await openDialog({
          title: `Copy ${currProjectName} to a new project`,
          maxWidth: 'md',
          yesText: 'Copy',
          initialValues: { newName: defaultName, newDir: defaultDir },
          render: ({ values, setValues }) => (
            <Stack spacing={2} sx={{ width: 560, maxWidth: '100%' }}>
              <TextField
                autoFocus
                fullWidth
                size='small'
                label='New project name'
                value={values.newName}
                onChange={(e) => setValues({ ...values, newName: e.target.value })}
                onClick={(e) => e.stopPropagation()}
              />
              <TextField
                fullWidth
                size='small'
                label='New files directory (leave empty for default)'
                value={values.newDir}
                onChange={(e) => setValues({ ...values, newDir: e.target.value })}
                onClick={(e) => e.stopPropagation()}
              />
            </Stack>
          ),
        });
        const newName = values?.newName.trim();
        if (confirmed && newName) {
          await copyProject(newName, values!.newDir.trim());
        }
      }}
    >
      <Stack direction={'row'} alignItems={'center'} spacing={0.5}>
        <ContentCopy fontSize="small" />
      </Stack>
    </ButtonTooltip>
    {DialogComponent}
  </>)
}
