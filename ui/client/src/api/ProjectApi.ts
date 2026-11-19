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

  // path starts at the project's files directory; "" is the directory itself.
  async files(projectName: string, path: string): Promise<any | null> {
    return await this.getJson('/project/files', { projectName, path }, { label: 'project files' });
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

  // The zip of the project's records; data files on disk are not in it.
  async export(projectName: string): Promise<ArrayBuffer | null> {
    return await this.postBytes('/project/export', {
      projectName,
    }, { label: `export project ${projectName}` });
  }

  // Loads an exported zip under the name it was exported with; answers the list and the project.
  async load(zipBytes: ArrayBuffer, fileName: string): Promise<any | null> {
    return await this.postBytesForJson('/project/load', zipBytes, { label: `load project ${fileName}` });
  }
}
