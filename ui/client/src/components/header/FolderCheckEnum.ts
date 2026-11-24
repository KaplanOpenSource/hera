// What the check of a project's files folder found so far.
export enum FolderCheckEnum {
  NoFolder = 'noFolder', // The project has no files folder set.
  Missing = 'missing', // A folder is set, but it is not on disk.
  Checking = 'checking', // The folder is on disk; we are still listing its files.
  Empty = 'empty', // The folder is on disk and holds nothing, or only caseConfiguration.json.
  NotEmpty = 'notEmpty', // The folder holds other files, or the check failed.
}
