// URL path of a project, optionally pointing at one of its documents.
export const projectPath = (projectName: string, docOid?: string): string => {
  const base = '/' + encodeURIComponent(projectName);
  if (docOid) {
    return `${base}/${docOid}`;
  }
  return base;
};
