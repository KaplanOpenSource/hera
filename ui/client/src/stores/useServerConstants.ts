import { create } from 'zustand';

interface ServerConstantStore {
  dataTypes: { [key: string]: string },
  setDataTypes: (dataTypes: { [key: string]: string }) => void,
}

// The data formats the server knows, read once at startup (see fetchInit).
export const useServerConstants = create<ServerConstantStore>((set) => ({
  dataTypes: {},
  setDataTypes: (dataTypes) => {
    set({ dataTypes })
  },
}));
