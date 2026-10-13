// A parameter's hermes name, written for people to read: "ProjectName" and
// "project_name" both become "Project Name".
export const friendlyParamName = (name: string): string => {
  const spaced = name
    .replace(/[_-]+/g, ' ')
    // Split "fooBar" into words, and "HTMLPage" into "HTML Page".
    .replace(/([a-z0-9])([A-Z])/g, '$1 $2')
    .replace(/([A-Z]+)([A-Z][a-z])/g, '$1 $2')
    .trim();
  const words: string[] = [];
  for (const word of spaced.split(/\s+/)) {
    if (word.length > 0) {
      words.push(word.charAt(0).toUpperCase() + word.slice(1));
    }
  }
  if (words.length === 0) {
    return name;
  }
  return words.join(' ');
};
