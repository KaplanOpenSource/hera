import { ProjectDocument, ProjectEntire } from "../shared/types";

export class ProjectObj {
  constructor(
    public data: ProjectEntire,
  ) {

  }

  documents(withoutConfig = true): DocumentObj[] {
    let docs = this.data.documents.map(d => new DocumentObj(d, this));
    if (withoutConfig) {
      docs = docs.filter(d => !d.isConfig());
    }
    return docs;
  }

  configDocument(): DocumentObj | undefined {
    return this.documents().filter(d => d.isConfig)[0];
  }
}

export class DocumentObj {
  constructor(
    public data: ProjectDocument,
    public project: ProjectObj,
  ) {

  }

  isConfig(): boolean {
    return this.data.type === this.project.data.name + '__config__';
  }
}