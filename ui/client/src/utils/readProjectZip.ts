import { unzipSync } from 'fflate';
import { Parser } from 'pickleparser';

export interface ProjectZipInfo {
  name: string,
  filesDirectory: string,
}

interface ZipRecord {
  projectName?: string,
  type?: string,
  desc?: { [key: string]: unknown },
}

// Each zip entry holds a list of records, or one record if exported without the progress bar.
const recordsOf = (bytes: Uint8Array): ZipRecord[] => {
  const parsed = new Parser().parse<ZipRecord | ZipRecord[]>(bytes);
  if (Array.isArray(parsed)) {
    return parsed;
  }
  return [parsed];
}

// The entry names in zip order, without unpacking any of them.
const entryNamesOf = (zipBytes: Uint8Array): string[] => {
  const names: string[] = [];
  unzipSync(zipBytes, {
    filter: (entry) => {
      names.push(entry.name);
      return false;
    },
  });
  return names;
}

// Reads the project name and files folder from a zip made by hera's Project.export.
export const readProjectZip = (zipBytes: Uint8Array): ProjectZipInfo => {
  const info: ProjectZipInfo = { name: '', filesDirectory: '' };
  // One entry at a time, stopping at the settings record, which is usually in the first.
  for (const entryName of entryNamesOf(zipBytes)) {
    const entryBytes = unzipSync(zipBytes, { filter: (entry) => entry.name === entryName })[entryName];
    for (const record of recordsOf(entryBytes)) {
      if (!info.name && record.projectName) {
        info.name = record.projectName;
      }
      if (record.projectName && record.type === `${record.projectName}__config__`) {
        info.filesDirectory = String(record.desc?.filesDirectory ?? '');
        return info;
      }
    }
  }
  return info;
}
