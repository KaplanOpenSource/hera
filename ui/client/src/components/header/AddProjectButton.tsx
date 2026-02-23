import { Add } from "@mui/icons-material";
import {
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogContentText,
  DialogTitle,
  TextField
} from "@mui/material";
import { ProjectEntire, ProjectName } from "@shared/types";
import { useRef, useState } from "react";
import { BooleanProperty } from "../../elements/BooleanProperty";
import { ButtonTooltip } from "../../elements/ButtonTooltip";
import { execPython } from "../../io/execPython";
import { useProjectStore } from "../../stores/useProjectStore";

export const AddProjectButton = ({ }) => {
  const [open, setOpen] = useState(false);
  const [name, setName] = useState('');
  const [filesDirectory, setFilesDirectory] = useState('');
  const [loadRepositories, setLoadRepositories] = useState(true);
  const { selectProject, setProjectNames, setCurrentProject } = useProjectStore();
  const inputRef = useRef();

  const doAddProject = async () => {
    let dirStr = `os.path.join(os.getcwd(), 'projects', '${name}')`;
    if (filesDirectory !== '') {
      dirStr = `'${filesDirectory}'`;
    }

    const { problem, data } = await execPython(`
import os
import json
from types import SimpleNamespace

from hera.utils.data.CLI import project_create
from hera.datalayer.project import getProjectList
from hera.datalayer import All
from hera import toolkitHome

project_create(SimpleNamespace(
  projectName='${name}',
  directory=${dirStr},
  loadRepositories=${loadRepositories ? 'True' : 'False'},
  overwrite=False))

projectNames = [{"name": proj} for proj in getProjectList()]

docs = All.getDocumentsAsDict('${name}', with_id=True)
project = {"name": '${name}', "documents": docs['documents']}

result = {"projectNames": projectNames, "project": project}
`)
    // table = toolkitHome.getToolkitTable('${name}')
    if (problem) {
      return;
    }

    // hera doesnt update new project immediately so we build names in front
    // await fetchProjectsNames();
    const details = (data?.project) as ProjectEntire;
    const names = (data?.projectNames || []) as ProjectName[];
    if (!names.find(x => x.name === name)) {
      names.push({ name });
    }
    setProjectNames(names);
    setCurrentProject(details);
    selectProject(name);
    setOpen(false);
  }

  return (<>
    <ButtonTooltip
      title='Add project'
      onClick={() => {
        setOpen(true);
        setName('');
        setTimeout(() => (inputRef.current as any)?.focus(), 0)
      }}
    >
      <Add />
    </ButtonTooltip>
    <Dialog
      open={open}
      onClose={() => setOpen(false)}
      onKeyDown={e => { if (e.code === 'Enter') doAddProject() }}
    >
      <DialogTitle>New Project</DialogTitle>
      <DialogContent>
        <DialogContentText>
          Adding a new project, please fill initial details
        </DialogContentText>
        <TextField
          inputRef={inputRef}
          // autoFocus
          required
          margin="dense"
          size="small"
          label="Project Name"
          fullWidth
          variant="outlined"
          value={name}
          onChange={(e) => setName(e.target.value)}
        />
        <TextField
          margin="dense"
          size="small"
          label={"Project Files Directory (leave empty for default)"}
          fullWidth
          variant="outlined"
          value={filesDirectory}
          onChange={(e) => setFilesDirectory(e.target.value)}
        />
        <BooleanProperty
          label="Load Repositories"
          value={loadRepositories}
          setValue={v => setLoadRepositories(v)}
        />
      </DialogContent>
      <DialogActions>
        <Button onClick={() => setOpen(false)}>
          Cancel
        </Button>
        <Button onClick={() => doAddProject()}>
          Add Project
        </Button>
      </DialogActions>
    </Dialog>
  </>)
}