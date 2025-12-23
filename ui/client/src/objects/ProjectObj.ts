import { ProjectDocument, ProjectEntire } from "../shared/types";

export class ProjectObj {
  constructor(
    public data: ProjectEntire,
  ) {
  }

  public get name(): string {
    return this.data.name;
  }

  public get allDocuments(): DocumentObj[] {
    return this.data.documents.map(d => new DocumentObj(d, this));
  }

  public get documents(): DocumentObj[] {
    return this.allDocuments.filter(d => !d.isConfig);
  }

  public get configDocument(): DocumentObj | undefined {
    return this.allDocuments.filter(d => d.isConfig)[0];
  }
}

export class DocumentObj {
  constructor(
    public data: ProjectDocument,
    public project: ProjectObj,
  ) {
  }

  public get docid(): string {
    return this.data._id.$oid;
  }

  public get isConfig(): boolean {
    return this.data.type === this.project.data.name + '__config__';
  }

  public get name(): string {
    return this.data.desc?.datasourceName || this.data.type || this.data._cls;

  }
}