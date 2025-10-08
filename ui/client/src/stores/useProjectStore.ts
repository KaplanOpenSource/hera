import { Project, ProjectName } from '@shared/types';
import { create } from 'zustand';

export const NO_PROJECT = "* NONE *";

interface ProjectStore {
  projectNames: ProjectName[]; // List of project names
  projectId: string;
  currProject: Project | null; // Current project
  setProjectNames: (names: ProjectName[]) => void; // Sets project names
  selectProject: (newProjectId: string) => void;
  setCurrentProject: (project: Project | null) => void; // Sets current project
}

export const useProjectStore = create<ProjectStore>((set) => ({
  projectNames: [],
  projectId: NO_PROJECT,
  currProject: null,
  setProjectNames: (names) => set({ projectNames: names }),
  selectProject: (newProjectId: string) => set({ projectId: newProjectId }),
  setCurrentProject: (project) => set({ currProject: project }),
}));
