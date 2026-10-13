import { Toolkit } from '@shared/types';
import { create } from 'zustand';
import { ToolkitObj } from '../objects/ToolkitObj';
import { useProjectStore } from './useProjectStore';

interface ToolkitStore {
  toolkits: Toolkit[];
  setToolkits: (val: Toolkit[]) => void;
  // The toolkit names actually used by the current project's documents.
  getProjectToolkitKeys: () => string[];
}

export const useToolkitStore = create<ToolkitStore>((set, get) => ({
  toolkits: [],
  setToolkits: (val) => {
    set({ toolkits: val })
  },
  getProjectToolkitKeys: () => {
    const { toolkits } = get();
    const documents = useProjectStore.getState().getProject()?.documents ?? [];
    const docToolkitNames = [...new Set(
      documents.map(d => d.toolkit).filter(Boolean) as string[]
    )];
    return toolkits
      .filter(t => docToolkitNames.some(dt =>
        new ToolkitObj(t).matches(dt)
      ))
      .map(t => t.toolkit);
  },
}));
