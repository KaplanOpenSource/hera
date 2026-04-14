import { ProjectDocument } from "../shared/types";
import { ProjectObj } from "./ProjectObj";

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

  public get isNotebook(): boolean {
    return this.data.type === 'notebook';
  }

  public get name(): string {
    return this.data.desc?.datasourceName || this.data.type || this.data._cls;
  }

  public get toolkit(): string | undefined {
    return this.data.desc.toolkit;
  }

  public get extDesc(): Record<string, any> {
    return { ...(this.data.desc || {}), type: this.data.type }
  }
}
