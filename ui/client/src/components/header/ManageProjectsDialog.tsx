import {
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  List,
  TextField,
  ThemeProvider,
  Typography
} from "@mui/material";
import { useEffect, useState } from "react";
import { projectApi } from "../../api/api";
import { ProjectFolder } from "../../shared/types";
import { DEFAULT_PROJECT } from "../../stores/useProjectStore";
import { useAppTheme } from "../../theme";
import { FolderCheckEnum } from "./FolderCheckEnum";
import { isFilesFolderEmpty } from "./isFilesFolderEmpty";
import { ManageProjectsRow } from "./ManageProjectsRow";

// Lists every project with its files folder, and notes the empty folders.
export const ManageProjectsDialog = ({
  onClose,
}: {
  onClose: () => void,
}) => {
  // Follow the app theme so the dialog isn't tinted by the dark header it opens from.
  const dialogTheme = useAppTheme();
  const [projects, setProjects] = useState<ProjectFolder[]>([]);
  const [folderChecks, setFolderChecks] = useState<{ [name: string]: FolderCheckEnum }>({});
  const [search, setSearch] = useState('');

  const setCheck = (name: string, check: FolderCheckEnum) => {
    setFolderChecks((prev) => {
      return { ...prev, [name]: check };
    });
  }

  const checkFolder = (name: string) => {
    isFilesFolderEmpty(name)
      .then((empty) => {
        if (empty) {
          setCheck(name, FolderCheckEnum.Empty);
        } else {
          setCheck(name, FolderCheckEnum.NotEmpty);
        }
      })
      // A failed check shows no note.
      .catch((error) => {
        console.error(`could not check the files folder of ${name}`, error);
        setCheck(name, FolderCheckEnum.NotEmpty);
      });
  }

  useEffect(() => {
    projectApi.folders().then((data) => {
      const shown: ProjectFolder[] = [];
      for (const project of (data?.projects ?? []) as ProjectFolder[]) {
        if (project.name !== DEFAULT_PROJECT) {
          shown.push(project);
        }
      }
      setProjects(shown);
      for (const project of shown) {
        if (!project.filesDirectory) {
          setCheck(project.name, FolderCheckEnum.NoFolder);
        } else if (!project.folderExists) {
          setCheck(project.name, FolderCheckEnum.Missing);
        } else {
          checkFolder(project.name);
        }
      }
    });
  }, []);

  const searchLower = search.trim().toLowerCase();
  const filtered: ProjectFolder[] = [];
  for (const project of projects) {
    if (project.name.toLowerCase().includes(searchLower)) {
      filtered.push(project);
    }
  }

  let countText = `${projects.length} projects`;
  if (projects.length === 1) {
    countText = '1 project';
  }

  return (
    <ThemeProvider theme={dialogTheme}>
      <Dialog
        open
        onClose={onClose}
        fullWidth
        maxWidth="sm"
      >
        <DialogTitle>Manage projects</DialogTitle>
        <DialogContent>
          <TextField
            label="Search"
            value={search}
            onChange={(e) => setSearch(e.target.value)}
            fullWidth
            margin="dense"
            size="small"
            slotProps={{ htmlInput: { 'data-testid': 'manage-projects-search' } }}
          />
          <Typography variant="body2" color="text.secondary" data-testid="manage-projects-count">
            {countText}
          </Typography>
          <List dense>
            {filtered.map((project) => (
              <ManageProjectsRow
                key={project.name}
                project={project}
                folderCheck={folderChecks[project.name] ?? FolderCheckEnum.Checking}
              />
            ))}
          </List>
        </DialogContent>
        <DialogActions>
          <Button onClick={onClose}>
            Close
          </Button>
        </DialogActions>
      </Dialog>
    </ThemeProvider>
  );
}
