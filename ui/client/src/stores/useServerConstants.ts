import { ProjectEntire, ProjectName, Toolkit } from '@shared/types';
import { useEffect } from 'react';
import { fetchPython } from '../io/fetchPython';
import { create } from 'zustand';

interface ServerConstantStore {
  dataTypes: { [key: string]: string },
  readAllConstants: () => Promise<void>,
}

export const useServerConstants = create<ServerConstantStore>((set) => ({
  dataTypes: {},
  readAllConstants: async () => {
    const { data } = await fetchPython({
      results: ['datatypes'],
      code: `
from hera import datalayer
datatypes = {key:value for key, value in vars(datalayer.datatypes).items()
    if not callable(value) and not key.startswith("__") and key[0] == key[0].upper()}
`,
    })
    if (data) {
      console.log('constants:', data)
      set({ dataTypes: data.datatypes || {} })
    }
  },
}));

export const ServerConstantReader = () => {
  useEffect(() => {
    useServerConstants.getState().readAllConstants();
  }, [])

  return null;
}