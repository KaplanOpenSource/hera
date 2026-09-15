import { NO_PROJECT, useProjectStore } from '../../../stores/useProjectStore';
import { WorkflowDesc } from '../../types';
import { WorkflowMutatorBase } from '../WorkflowMutatorBase';

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

// Phase: seed every empty project-name field in the document's desc with the
// current project. Covers a field the user just renamed to ProjectName, one
// seeded by adding a node, and one that was already there and empty.
export class FillProjectNameMutator extends WorkflowMutatorBase {
  readonly name = 'fillProjectName';

  mutate(desc: WorkflowDesc): WorkflowDesc {
    const projectName = useProjectStore.getState().currProjectName;
    if (!projectName || projectName === NO_PROJECT) {
      return desc;
    }
    return fillProjectName(desc, projectName);
  }
}
