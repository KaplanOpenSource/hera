import { Project, ProjectName } from '@shared/types';
import { create } from 'zustand';

interface ProjectStore {
  projectNames: ProjectName[]; // List of project names
  currProject: Project | null; // Current project
  setProjectNames: (names: ProjectName[]) => void; // Sets project names
  setCurrentProject: (project: Project | null) => void; // Sets current project
}

export const useProjectStore = create<ProjectStore>((set) => ({
  projectNames: [],
  currProject: null,
  setProjectNames: (names) => set({ projectNames: names }),
  setCurrentProject: (project) => set({ currProject: project }),
}));
