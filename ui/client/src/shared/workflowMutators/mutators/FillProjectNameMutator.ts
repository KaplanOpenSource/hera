// True for a field name that names a project: `ProjectName`, `projectName`,
// `projectname` - any casing. Hera/Hermes nodes use `ProjectName`.
export const isProjectNameKey = (key: string): boolean => {
  return /^projectname$/i.test(key);
};

// Empty means no value yet: unset, null, or the empty string.
const isEmpty = (value: any): boolean => {
  return value === undefined || value === null || value === '';
};

// Walks a value and fills every empty project-name field at any depth, in
// objects and in arrays alike. A field that already holds a value is left
// alone, so a user can point it at another project. Never adds a field.
export const fillProjectName = (value: any, projectName: string): any => {
  if (Array.isArray(value)) {
    return value.map((item) => {
      return fillProjectName(item, projectName);
    });
  }
  if (value === null || typeof value !== 'object') {
    return value;
  }
  const next: { [key: string]: any } = {};
  for (const [key, item] of Object.entries(value)) {
    if (isProjectNameKey(key) && isEmpty(item)) {
      next[key] = projectName;
    } else {
      next[key] = fillProjectName(item, projectName);
    }
  }
  return next;
};
