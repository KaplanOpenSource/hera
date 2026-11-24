import {
  Button,
  Dialog,
  DialogActions,
  DialogContent,
  DialogTitle,
  Table,
  TableBody,
  TableCell,
  TableHead,
  TableRow,
  TextField,
  ThemeProvider,
  Typography
} from "@mui/material";
import { useEffect, useState } from "react";
import { projectApi } from "../../api/api";
import { ProjectFolder } from "../../shared/types";
import { DEFAULT_PROJECT } from "../../stores/useProjectStore";
import { useAppTheme } from "../../theme";
import { ExportProgress } from "./ExportProgress";
import { ExportProjectsButton } from "./ExportProjectsButton";
import { ExportStep } from "./ExportStep";
import { FolderCheckEnum } from "./FolderCheckEnum";
import { isFilesFolderEmpty } from "./isFilesFolderEmpty";
import { ManageProjectsRow } from "./ManageProjectsRow";
import { MarkAllCheckbox } from "./MarkAllCheckbox";

// Lists every project with its files folder, notes the empty folders, and exports the marked ones.
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
  const [marks, setMarks] = useState<{ [name: string]: boolean }>({});
  const [exportStep, setExportStep] = useState<ExportStep | null>(null);
  const [lastClicked, setLastClicked] = useState<string | null>(null);

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

  // The shown names from the last clicked row to this one, both included.
  const namesFromLastClicked = (name: string) => {
    const shownNames: string[] = [];
    for (const project of filtered) {
      shownNames.push(project.name);
    }
    const from = shownNames.indexOf(lastClicked ?? '');
    const to = shownNames.indexOf(name);
    if (from < 0) {
      return [name];
    }
    return shownNames.slice(Math.min(from, to), Math.max(from, to) + 1);
  }

  // With shift, every row between the last clicked one and this one takes this row's new mark.
  const toggleMark = (name: string, shift: boolean) => {
    const mark = !marks[name];
    let names = [name];
    if (shift) {
      names = namesFromLastClicked(name);
    }
    setMarks((prev) => {
      const next = { ...prev };
      for (const each of names) {
        next[each] = mark;
      }
      return next;
    });
    setLastClicked(name);
  }

  let countText = `${projects.length} projects`;
  if (projects.length === 1) {
    countText = '1 project';
  }

  const markedNames: string[] = [];
  for (const project of projects) {
    if (marks[project.name]) {
      markedNames.push(project.name);
    }
  }

  const markAll = (mark: boolean) => {
    const next: { [name: string]: boolean } = {};
    if (mark) {
      for (const project of projects) {
        next[project.name] = true;
      }
    }
    setMarks(next);
  }

  return (
    <ThemeProvider theme={dialogTheme}>
      <Dialog
        open
        onClose={onClose}
        fullWidth
        maxWidth="lg"
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
          <Typography variant="body2" color="text.secondary" data-testid="manage-projects-marked">
            {`${markedNames.length} of ${projects.length} marked`}
          </Typography>
          {/* Fixed layout so long folder paths wrap instead of widening the dialog. */}
          <Table size="small" sx={{ tableLayout: 'fixed' }}>
            <TableHead>
              <TableRow>
                <TableCell padding="checkbox">
                  <MarkAllCheckbox
                    markedCount={markedNames.length}
                    total={projects.length}
                    onMarkAll={markAll}
                  />
                </TableCell>
                <TableCell sx={{ width: '30%' }}>Name</TableCell>
                <TableCell sx={{ width: '45%' }}>Files folder</TableCell>
                <TableCell sx={{ width: '25%' }}>Note</TableCell>
              </TableRow>
            </TableHead>
            <TableBody>
              {filtered.map((project) => (
                <ManageProjectsRow
                  key={project.name}
                  project={project}
                  folderCheck={folderChecks[project.name] ?? FolderCheckEnum.Checking}
                  marked={!!marks[project.name]}
                  onToggle={(shift) => toggleMark(project.name, shift)}
                />
              ))}
            </TableBody>
          </Table>
        </DialogContent>
        {exportStep && <ExportProgress step={exportStep} />}
        <DialogActions>
          <ExportProjectsButton names={markedNames} onProgress={setExportStep} />
          <Button onClick={onClose}>
            Close
          </Button>
        </DialogActions>
      </Dialog>
    </ThemeProvider>
  );
}
