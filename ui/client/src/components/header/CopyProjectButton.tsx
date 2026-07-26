import { ContentCopy } from "@mui/icons-material";
import { Box, Checkbox, FormControlLabel, Stack, TextField, Typography } from "@mui/material";
import { useNavigate } from "react-router-dom";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { useDialog } from "../../elements/useDialog";
import { fetchPython } from "../../io/fetchPython";
import { ProjectName } from "../../shared/types";
import { useProjectStore } from "../../stores/useProjectStore";

interface CopyProjectValues {
  newName: string;
  newDir: string;
  selectedFolders: string[];
}

export const CopyProjectButton = ({
  onCopied,
}: {
  onCopied?: () => void,
}) => {
  const { openDialog, DialogComponent } = useDialog<CopyProjectValues>();
  const { currProjectName, setProjectNames } = useProjectStore();
  const navigate = useNavigate();

  // Default the copy into the source's "projects" folder (the nearest 'projects' ancestor,
  // or a projects/ subfolder if the source isn't under one), named after the copy.
  const srcFilesDir = useProjectStore.getState().getProject()?.configDocument?.data.desc.filesDirectory ?? '';
  const defaultName = `${currProjectName}_copy`;
  const cleanDir = srcFilesDir.replace(/\/+$/, '');
  const parts = cleanDir.split('/');
  const projectsIdx = parts.lastIndexOf('projects');
  const projectsDir = projectsIdx >= 0
    ? parts.slice(0, projectsIdx + 1).join('/')
    : (cleanDir ? `${cleanDir}/projects` : '');
  const defaultDir = projectsDir ? `${projectsDir}/${defaultName}` : '';

  // File resources of the current project's file-backed documents.
  const sourceFiles = (useProjectStore.getState().getProject()?.documents ?? [])
    .map(d => d.data.resource)
    .filter((r): r is string => typeof r === 'string' && r.includes('/'));

  // Distinct folders those files live in.
  const sourceFolders = [...new Set(
    sourceFiles.map(r => r.slice(0, r.lastIndexOf('/')))
  )].sort();

  // Only folders inside the project's files directory can be copied; the rest are collapsed into one "outside" row.
  const underProject = (p: string) => !!srcFilesDir && (p === srcFilesDir || p.startsWith(srcFilesDir + '/'));
  const inProjectFolders = sourceFolders.filter(underProject);
  const hasOutside = sourceFolders.some(f => !underProject(f));

  const copyProject = async (newName: string, newDir: string, selectedFolders: string[]) => {
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
selectedFolders = ${JSON.stringify(selectedFolders)}
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
    # Copy single files whose folder was selected (all selected folders live inside the project);
    # everything else keeps its original resource path.
    if isinstance(res, str) and res:
        resAbs = os.path.abspath(res)
        if os.path.isfile(resAbs) and os.path.dirname(resAbs) in selectedFolders:
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
        console.log('Project files:', sourceFiles);
        const { confirmed, values } = await openDialog({
          title: `Copy ${currProjectName} to a new project`,
          maxWidth: 'md',
          yesText: 'Copy',
          initialValues: { newName: defaultName, newDir: defaultDir, selectedFolders: inProjectFolders },
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
              {inProjectFolders.length === 0 && !hasOutside
                ? <Typography variant='caption' color='text.secondary'>No file-backed documents.</Typography>
                : (
                  <Box>
                    <Typography variant='caption' color='text.secondary' sx={{ display: 'block', mb: 0.5 }}>
                      Copy files from folders:
                    </Typography>
                    <Box sx={{ maxHeight: 200, overflowY: 'auto' }}>
                      {inProjectFolders.map(f => (
                        <FormControlLabel
                          key={f}
                          sx={{ display: 'flex', m: 0 }}
                          control={
                            <Checkbox
                              size='small'
                              sx={{ py: 0.25 }}
                              checked={values.selectedFolders.includes(f)}
                              onChange={(e) => setValues({
                                ...values,
                                selectedFolders: e.target.checked
                                  ? [...values.selectedFolders, f]
                                  : values.selectedFolders.filter(x => x !== f),
                              })}
                            />
                          }
                          label={<Box component='span' sx={{ fontFamily: 'monospace', fontSize: 13 }}>{f}</Box>}
                        />
                      ))}
                    </Box>
                    {hasOutside && (
                      <Typography variant='caption' color='text.secondary' sx={{ mt: 0.5, display: 'block' }}>
                        Some files are outside the project folder and will not be copied.
                      </Typography>
                    )}
                  </Box>
                )
              }
            </Stack>
          ),
        });
        const newName = values?.newName.trim();
        if (confirmed && newName) {
          await copyProject(newName, values!.newDir.trim(), values!.selectedFolders);
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
