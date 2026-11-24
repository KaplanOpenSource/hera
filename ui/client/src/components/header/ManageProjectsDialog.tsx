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
import { DeleteProjectsButton } from "./DeleteProjectsButton";
import { ExportProjectsButton } from "./ExportProjectsButton";
import { EmptyProjectsToggle } from "./EmptyProjectsToggle";
import { isEmptyProject } from "./isEmptyProject";
import { ManageProjectsRow } from "./ManageProjectsRow";
import { MarkAllCheckbox } from "./MarkAllCheckbox";

// Lists every project with its files folder, notes the folder state, and exports or deletes the marked ones.
export const ManageProjectsDialog = ({
  onClose,
}: {
  onClose: () => void,
}) => {
  // Follow the app theme so the dialog isn't tinted by the dark header it opens from.
  const dialogTheme = useAppTheme();
  const [projects, setProjects] = useState<ProjectFolder[]>([]);
  const [search, setSearch] = useState('');
  const [onlyEmpty, setOnlyEmpty] = useState(false);
  const [marks, setMarks] = useState<{ [name: string]: boolean }>({});
  const [lastClicked, setLastClicked] = useState<string | null>(null);

  useEffect(() => {
    projectApi.folders().then((data) => {
      const shown: ProjectFolder[] = [];
      for (const project of (data?.projects ?? []) as ProjectFolder[]) {
        if (project.name !== DEFAULT_PROJECT) {
          shown.push(project);
        }
      }
      setProjects(shown);
    });
  }, []);

  const searchLower = search.trim().toLowerCase();
  const filtered: ProjectFolder[] = [];
  for (const project of projects) {
    const matchesSearch = project.name.toLowerCase().includes(searchLower);
    const matchesEmpty = !onlyEmpty || isEmptyProject(project);
    if (matchesSearch && matchesEmpty) {
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
  const markedProjects: ProjectFolder[] = [];
  for (const project of projects) {
    if (marks[project.name]) {
      markedNames.push(project.name);
      markedProjects.push(project);
    }
  }

  // Drops the deleted projects' rows and marks.
  const removeDeleted = (deleted: string[]) => {
    setProjects((prev) => {
      const left: ProjectFolder[] = [];
      for (const project of prev) {
        if (!deleted.includes(project.name)) {
          left.push(project);
        }
      }
      return left;
    });
    setMarks((prev) => {
      const next = { ...prev };
      for (const name of deleted) {
        delete next[name];
      }
      return next;
    });
    if (lastClicked && deleted.includes(lastClicked)) {
      setLastClicked(null);
    }
  }

  let shownMarkedCount = 0;
  for (const project of filtered) {
    if (marks[project.name]) {
      shownMarkedCount += 1;
    }
  }

  // Marking leaves exactly the shown rows marked, so a hidden row is never marked by it; clearing clears all.
  const markAll = (mark: boolean) => {
    const next: { [name: string]: boolean } = {};
    if (mark) {
      for (const project of filtered) {
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
          <EmptyProjectsToggle
            onlyEmpty={onlyEmpty}
            onChange={setOnlyEmpty}
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
                    markedCount={shownMarkedCount}
                    total={filtered.length}
                    onMarkAll={markAll}
                  />
                </TableCell>
                <TableCell sx={{ width: '22%' }}>Name</TableCell>
                <TableCell sx={{ width: '35%' }}>Files folder</TableCell>
                <TableCell sx={{ width: '10%' }}>Documents</TableCell>
                <TableCell sx={{ width: '18%' }}>Note</TableCell>
                <TableCell sx={{ width: '15%' }} />
              </TableRow>
            </TableHead>
            <TableBody>
              {filtered.map((project) => (
                <ManageProjectsRow
                  key={project.name}
                  project={project}
                  marked={!!marks[project.name]}
                  onToggle={(shift) => toggleMark(project.name, shift)}
                  onSwitched={onClose}
                />
              ))}
            </TableBody>
          </Table>
        </DialogContent>
        <DialogActions>
          <ExportProjectsButton
            names={markedNames}
          />
          <DeleteProjectsButton
            projects={markedProjects}
            onDeleted={removeDeleted}
          />
          <Button onClick={onClose}>
            Close
          </Button>
        </DialogActions>
      </Dialog>
    </ThemeProvider>
  );
}
