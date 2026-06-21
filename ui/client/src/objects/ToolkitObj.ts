import { Toolkit } from "../shared/types";

export class ToolkitObj {
  public readonly registered: boolean;

  constructor(
    public data: Toolkit,
    registered = true,
  ) {
    this.registered = registered;
  }

  static unregistered(toolkitName: string): ToolkitObj {
    return new ToolkitObj({ toolkit: toolkitName, cls: '' }, false);
  }

  get toolkit(): string {
    return this.data.toolkit;
  }

  get cls(): string {
    return this.data.cls;
  }

  get className(): string {
    return this.data.cls.split('.').pop() ?? '';
  }

  get shortName(): string {
    return this.className.replace(/Toolkit$/i, '') || this.data.toolkit;
  }

  get description(): string | undefined {
    return this.data.description;
  }

  get source(): string | undefined {
    return this.data.source;
  }

  get type(): string | undefined {
    return this.data.type;
  }

  get repositoryName(): string | undefined {
    return this.data.repositoryName;
  }

  get version(): string | undefined {
    return this.data.version;
  }

  matches(name: string | undefined): boolean {
    if (!name) return !this.registered;
    return this.toolkit === name || this.shortName === name || this.className === name;
  }
}
