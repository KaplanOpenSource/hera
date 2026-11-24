// What a project's files folder holds, as the server's folders list says.
export enum FolderStateEnum {
  NoFolder = 'noFolder', // The project has no files folder set.
  Missing = 'missing', // A folder is set, but nothing is at that path.
  NotAFolder = 'notAFolder', // Something is at that path, but it is a file.
  Empty = 'empty', // The folder holds nothing.
  OnlyCaseConfig = 'onlyCaseConfig', // The folder holds only caseConfiguration.json.
  HasFiles = 'hasFiles', // The folder holds other files or subfolders.
  Unreadable = 'unreadable', // The folder is there but cannot be listed.
}
