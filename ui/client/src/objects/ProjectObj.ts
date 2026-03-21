import { ProjectEntire } from "../shared/types";
import { DocumentObj } from "./DocumentObj";

export { DocumentObj } from "./DocumentObj";

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

  public get documentIds(): Set<string> {
    return new Set(this.allDocuments.map(d => d.docid));
  }

  public get configDocument(): DocumentObj | undefined {
    return this.allDocuments.filter(d => d.isConfig)[0];
  }
}
