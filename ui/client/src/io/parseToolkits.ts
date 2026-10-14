import { Toolkit } from "@shared/types";

// Turns toolkit documents into the toolkits the UI shows.
export const parseToolkits = (toolkitDocs: any[]): Toolkit[] => {
  return toolkitDocs.map((d: any) => {
    const desc = d.desc;
    const t = {
      toolkit: d.toolkit,
      cls: desc.classpath,
    } as Toolkit;
    for (const field of [
      'description',
      'source',
      'type',
      'repositoryName',
      'version',
    ] as const) {
      if (desc[field]) t[field] = desc[field];
    }
    return t;
  });
};
