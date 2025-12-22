import { ProjectEntire, ProjectName, Toolkit } from '@shared/types';
import { useEffect } from 'react';
import { execPython } from '../io/execPython';
import { create } from 'zustand';

interface ServerConstantStore {
  dataTypes: { [key: string]: string },
  readAllConstants: () => Promise<void>,
}

export const useServerConstants = create<ServerConstantStore>((set) => ({
  dataTypes: {},
  readAllConstants: async () => {
    const { data } = await execPython(`
from hera import datalayer
datatypes = {key:value for key, value in vars(datalayer.datatypes).items()
    if not callable(value) and not key.startswith("__") and key[0] == key[0].upper()}
result = {'datatypes': datatypes}
`)
    if (data) {
      console.log('constants:', data)
      set({ dataTypes: data['datatypes'] || {} })
    }
  },
}));

export const ServerConstantReader = () => {
  console.log(1)
  useEffect(() => {
    useServerConstants.getState().readAllConstants();
  }, [])

  return null;
}