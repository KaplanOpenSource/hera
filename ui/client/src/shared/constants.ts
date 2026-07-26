
export const FORBIDDEN_FIELDS = ['_id', '_cls', 'projectName'];

// The document's dataFormat field, which gets its own editor instead of the
// generic value editor / type chip.
export const DATA_FORMAT_FIELD = 'dataFormat';

// The document's desc field (holds the descriptive sub-fields).
export const DESC_FIELD = 'desc';

// desc's filesDirectory sub-field, which is shown read-only (can't be renamed).
export const FILES_DIRECTORY_FIELD = 'filesDirectory';
