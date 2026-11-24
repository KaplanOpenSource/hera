// Where a multi-project export is: how many are done, out of how many, and which one runs now.
export interface ExportStep {
  done: number,
  total: number,
  current: string,
}
