import { SECTION_DATASOURCE, SECTION_EXPERIMENT } from "./RepoJsonMerger";

export const repoIsToolkit = (treePath: string[]) =>
  treePath.length === 1 && treePath[0].toLowerCase() !== SECTION_EXPERIMENT;

export const repoIsExperimentSection = (treePath: string[]) =>
  treePath.length === 1 && treePath[0].toLowerCase() === SECTION_EXPERIMENT;

export const repoIsExperiment = (treePath: string[]) =>
  treePath.length === 2 && treePath[0].toLowerCase() === SECTION_EXPERIMENT;

export const repoIsDatasource = (treePath: string[]) =>
  treePath.length === 3 && treePath[1].toLowerCase() === SECTION_DATASOURCE;

