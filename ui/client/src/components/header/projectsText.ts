// "1 project" or "3 projects".
export const projectsText = (count: number) => {
  if (count === 1) {
    return '1 project';
  }
  return `${count} projects`;
}
