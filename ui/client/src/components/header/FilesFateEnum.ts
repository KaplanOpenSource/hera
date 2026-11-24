// What a multi delete does to a project's files folder.
export enum FilesFateEnum {
  // The folder is empty or holds only the case config, so the empty folder is removed.
  Removed = 'removed',
  // The folder has files, is a file, or cannot be read, so it stays.
  Kept = 'kept',
  // No folder is set or it is missing, so there is nothing to delete.
  Nothing = 'nothing',
}
