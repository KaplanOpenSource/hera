import { ProjectDocument, ProjectEntire } from '../../../shared/types';

// The project's config document (type "<name>__config__"), if it has one.
export const findConfigDoc = (project: ProjectEntire | null): ProjectDocument | undefined => {
  return project?.documents.find(d => d.type === `${project.name}__config__`);
};
