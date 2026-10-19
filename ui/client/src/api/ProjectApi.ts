import { ApiCall } from './ApiCall';

// Every call to /project/*. Each method answers data, or null when the call
// failed.
export class ProjectApi extends ApiCall {
  async names(): Promise<any | null> {
    return await this.getJson('/project/names', {}, { label: 'project names' });
  }

  // label varies by caller; silent is for the periodic auto-reload, no toasts.
  async documents(projectName: string, label: string, silent = false): Promise<any | null> {
    return await this.getJson('/project/documents', { projectName }, { label, silent });
  }

  async toolkits(): Promise<any | null> {
    return await this.getJson('/project/toolkits', {}, { label: 'toolkits' });
  }

  async create(projectName: string, filesDirectory: string, loadRepositories: boolean): Promise<any | null> {
    return await this.postJson('/project/create', {
      projectName,
      filesDirectory,
      loadRepositories,
    }, { label: `create project ${projectName}` });
  }

  async detectNotebooks(projectName: string, filesDirectory: string): Promise<any | null> {
    return await this.postJson('/project/detectNotebooks', {
      projectName,
      filesDirectory,
    }, { label: 'detect notebooks' });
  }

  async delete(projectName: string, deleteFiles: boolean): Promise<any | null> {
    return await this.postJson('/project/delete', {
      projectName,
      deleteFiles,
    }, { label: `delete project ${projectName}` });
  }
}
