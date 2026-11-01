import { Description, Folder } from '@mui/icons-material';
import { Box, Stack, Typography } from '@mui/material';
import { SimpleTreeView, TreeItem } from '@mui/x-tree-view';
import { ReactNode, useEffect, useState } from 'react';
import { projectApi } from '../../api/api';
import { ProjectObj } from '../../objects/ProjectObj';

interface FileEntry {
  name: string,
  isFolder: boolean,
}

// Shown under a folder until its contents arrive, so it gets an expand arrow.
const LOADING_SUFFIX = '/__loading__';

// One folder's contents, by path from the project folder; "" is the folder itself.
type FolderContents = { [path: string]: FileEntry[] };

const fetchFolder = async (projectName: string, path: string): Promise<FileEntry[]> => {
  const data = await projectApi.files(projectName, path);
  return data?.files ?? [];
};

const entryLabel = (entry: FileEntry) => {
  let icon = <Description sx={{ fontSize: 18, color: 'text.secondary' }} />;
  if (entry.isFolder) {
    icon = <Folder sx={{ fontSize: 18 }} />;
  }
  return (
    <Stack direction='row' alignItems='center' spacing={1}>
      {icon}
      <Typography>{entry.name}</Typography>
    </Stack>
  );
};

// The project folder as a tree. A folder is read from the server only when opened.
export const ProjectFilesView = ({
  project,
}: {
  project: ProjectObj,
}) => {
  const filesDirectory: string = project.configDocument?.data.desc.filesDirectory ?? '';
  const [contents, setContents] = useState<FolderContents>({});

  // Start again when the project or its folder changes.
  useEffect(() => {
    setContents({});
    if (!filesDirectory) {
      return;
    }
    let current = true;
    fetchFolder(project.name, '').then(entries => {
      if (current) {
        setContents({ '': entries });
      }
    });
    return () => {
      current = false;
    };
  }, [project.name, filesDirectory]);

  const openFolder = async (path: string) => {
    if (contents[path]) {
      return;
    }
    const entries = await fetchFolder(project.name, path);
    setContents(prev => ({ ...prev, [path]: entries }));
  };

  const renderFolder = (path: string): ReactNode[] => {
    const items: ReactNode[] = [];
    for (const entry of contents[path] ?? []) {
      let childPath = entry.name;
      if (path) {
        childPath = `${path}/${entry.name}`;
      }
      let children: ReactNode = null;
      if (entry.isFolder && contents[childPath]) {
        children = renderFolder(childPath);
      } else if (entry.isFolder) {
        children = <TreeItem itemId={childPath + LOADING_SUFFIX} label='Loading...' disabled />;
      }
      items.push(
        <TreeItem key={childPath} itemId={childPath} label={entryLabel(entry)}>
          {children}
        </TreeItem>
      );
    }
    return items;
  };

  let body: ReactNode = (
    <SimpleTreeView
      key={project.name + filesDirectory}
      onItemExpansionToggle={(_e, itemId, isExpanded) => isExpanded && openFolder(itemId)}
    >
      {renderFolder('')}
    </SimpleTreeView>
  );
  if (!filesDirectory) {
    body = <Typography color='text.secondary'>This project has no files directory.</Typography>;
  } else if (contents[''] && contents[''].length === 0) {
    body = <Typography color='text.secondary'>The folder is empty or not found.</Typography>;
  }

  return (
    <Box sx={{ p: 2, height: '100%', overflow: 'auto' }}>
      <Typography variant='subtitle2' color='text.secondary' sx={{ mb: 1 }}>
        {filesDirectory}
      </Typography>
      {body}
    </Box>
  );
};
