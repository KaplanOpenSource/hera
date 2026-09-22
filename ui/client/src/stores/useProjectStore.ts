import { ProjectDocument, ProjectEntire } from '@shared/types';
import { create } from 'zustand';
import { ProjectObj } from '../objects/ProjectObj';

export const NO_PROJECT = "* NONE *";
export const DEFAULT_PROJECT = "defaultProject";
export const EMPTY_NAME_PROJECT = "* Empty Name *";

interface ProjectStore {
  currProjectName: string;
  currProject: ProjectEntire | null; // Current project
  editedDocs: { [docid: string]: ProjectDocument }; // Unsaved edits, by document id
  selectProject: (newProjectId: string) => void;
  setCurrentProject: (project: ProjectEntire | null) => void; // Sets current project
  setEditedDoc: (docid: string, doc: ProjectDocument | null) => void; // null drops the edits
  getProject: () => ProjectObj | null;
}

export const useProjectStore = create<ProjectStore>((set, get) => ({
  currProjectName: NO_PROJECT,
  currProject: null,
  editedDocs: {},
  selectProject: (newProjectName: string) => {
    set({ currProjectName: newProjectName })
  },
  setCurrentProject: (project) => {
    const prevEdited = get().editedDocs;
    if (Object.keys(prevEdited).length === 0) {
      set({ currProject: project })
      return;
    }
    // Drop edits for documents the new project doesn't have.
    const ids = project ? new ProjectObj(project).documentIds : new Set<string>();
    const editedDocs: { [docid: string]: ProjectDocument } = {};
    for (const [docid, doc] of Object.entries(prevEdited)) {
      if (ids.has(docid)) {
        editedDocs[docid] = doc;
      }
    }
    set({ currProject: project, editedDocs })
  },
  setEditedDoc: (docid, doc) => {
    const editedDocs = { ...get().editedDocs };
    if (doc) {
      editedDocs[docid] = doc;
    } else {
      delete editedDocs[docid];
    }
    set({ editedDocs })
  },
  getProject: () => {
    const { currProject } = get();
    return currProject ? new ProjectObj(currProject) : null;
  },
}));
