import { DocumentApi } from './DocumentApi';
import { NodeApi } from './NodeApi';
import { ProjectApi } from './ProjectApi';
import { Reporter } from './Reporter';
import { RepositoryApi } from './RepositoryApi';
import { ServerApi } from './ServerApi';
import { WorkflowApi } from './WorkflowApi';

export const projectApi = new ProjectApi();
export const documentApi = new DocumentApi();
export const repositoryApi = new RepositoryApi();
export const workflowApi = new WorkflowApi();
export const nodeApi = new NodeApi();
export const serverApi = new ServerApi();

const all = [projectApi, documentApi, repositoryApi, workflowApi, nodeApi, serverApi];

// Called once at startup. Until then the calls work and say nothing.
export const installReporter = (reporter: Reporter) => {
  for (const api of all) {
    api.setReporter(reporter);
  }
};
