export class RepoJsonMerger {
  merged: { [key: string]: any } = {};
  overrides: { [path: string]: string[] } = {};

  constructor(repoJsons: { [path: string]: { [key: string]: any } }) {
    for (const [sourcePath, doc] of Object.entries(repoJsons)) {
      for (const [toolkit, sections] of Object.entries(doc)) {
        if (!this.merged[toolkit]) {
          this.merged[toolkit] = {};
        }
        if (toolkit.toLowerCase() === 'experiment') {
          const dsSection = Object.entries(sections as { [key: string]: any })
            .find(([key]) => key.toLowerCase() === 'datasource');
          const dsName = dsSection
            && typeof dsSection[1] === 'object'
            && dsSection[1] !== null
            && Object.keys(dsSection[1]).length > 0
            ? Object.keys(dsSection[1])[0]
            : null;
          if (dsName) {
            if (!this.merged[toolkit][dsName]) {
              this.merged[toolkit][dsName] = {};
            }
            for (const [section, sectionData] of Object.entries(sections as { [key: string]: any })) {
              if (!this.merged[toolkit][dsName][section]) {
                this.merged[toolkit][dsName][section] = {};
              }
              this.mergeInto(
                this.merged[toolkit][dsName][section],
                sectionData as { [key: string]: any },
                sourcePath,
                `${toolkit}/${dsName}/${section}`,
              );
            }
            continue;
          }
        }
        for (const [section, sectionData] of Object.entries(sections as { [key: string]: any })) {
          if (!this.merged[toolkit][section]) {
            this.merged[toolkit][section] = {};
          }
          this.mergeInto(
            this.merged[toolkit][section],
            sectionData as { [key: string]: any },
            sourcePath,
            `${toolkit}/${section}`,
          );
        }
      }
    }
  }

  private mergeInto(
    target: { [key: string]: any },
    source: { [key: string]: any },
    sourcePath: string,
    pathPrefix: string,
  ) {
    for (const [key, val] of Object.entries(source)) {
      const fullPath = `${pathPrefix}/${key}`;
      if (key in target) {
        this.diffAndRecord(target[key], val, sourcePath, fullPath);
      }
      target[key] = val;
    }
  }

  private diffAndRecord(
    oldVal: any,
    newVal: any,
    sourcePath: string,
    path: string,
  ) {
    if (typeof oldVal === 'object' && oldVal !== null
      && typeof newVal === 'object' && newVal !== null
      && !Array.isArray(oldVal) && !Array.isArray(newVal)) {
      const allKeys = new Set([...Object.keys(oldVal), ...Object.keys(newVal)]);
      for (const key of allKeys) {
        if (!(key in oldVal) || !(key in newVal)) {
          this.recordOverride(`${path}/${key}`, sourcePath);
        } else {
          this.diffAndRecord(oldVal[key], newVal[key], sourcePath, `${path}/${key}`);
        }
      }
      return;
    }
    if (JSON.stringify(oldVal) !== JSON.stringify(newVal)) {
      this.recordOverride(path, sourcePath);
    }
  }

  private recordOverride(path: string, sourcePath: string) {
    if (this.overrides[path]) {
      this.overrides[path].push(sourcePath);
    } else {
      this.overrides[path] = [sourcePath];
    }
  }
}
