// The tree row a click landed on: the newly selected one, or the only one left selected.
export const clickedRow = (ids: string[], prevIds: string[]) => {
  const added = ids.filter(id => !prevIds.includes(id));
  if (added.length > 0) {
    return added[added.length - 1];
  }
  if (ids.length === 1) {
    return ids[0];
  }
  return undefined;
};
