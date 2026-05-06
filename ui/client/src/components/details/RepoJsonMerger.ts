export class RepoJsonMerger {
  merged: { [key: string]: any } = {};
  overrides = new Map<string, string[]>();

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
      if (key in target && JSON.stringify(target[key]) !== JSON.stringify(val)) {
        const existing = this.overrides.get(fullPath);
        if (existing) {
          existing.push(sourcePath);
        } else {
          this.overrides.set(fullPath, [sourcePath]);
        }
      }
      target[key] = val;
    }
  }
}
