// True for a field named ProjectName in any casing.
export const isProjectNameKey = (key: string): boolean => {
  return /^projectname$/i.test(key);
};

// Empty means no value yet: unset, null, or the empty string.
const isEmpty = (value: any): boolean => {
  return value === undefined || value === null || value === '';
};

// Fills every empty project-name field at any depth, in objects and arrays. A
// field that already holds a value is left alone, and none is ever added.
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
