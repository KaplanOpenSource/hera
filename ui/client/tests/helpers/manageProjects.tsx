import { render, screen, fireEvent, within, act } from '@testing-library/react';
import { MemoryRouter, useLocation } from 'react-router-dom';
import { ManageProjectsButton } from '../../src/components/header/ManageProjectsButton';

// Shared by the manage projects tests: the folders answer, opening the dialog, and finding rows.

export const FOLDERS = {
  projects: [
    { name: 'defaultProject', filesDirectory: '/data/default', folderState: 'hasFiles', documentCount: 1 },
    { name: 'Alpha', filesDirectory: '/data/alpha', folderState: 'hasFiles', documentCount: 3 },
    { name: 'Beta', filesDirectory: '/data/beta', folderState: 'hasFiles', documentCount: 0 },
    { name: 'Gamma', filesDirectory: '/data/gamma', folderState: 'hasFiles', documentCount: 12 },
  ],
};

// Shows the current route, so tests can see where a switch went.
const LocationShower = () => {
  const location = useLocation();
  return <span data-testid="location">{location.pathname}</span>;
};

export const openDialog = async () => {
  render(<MemoryRouter><ManageProjectsButton /><LocationShower /></MemoryRouter>);
  const btn = within(screen.getByLabelText('Manage projects')).getByRole('button');
  await act(async () => { fireEvent.click(btn); });
  await screen.findByText('Manage projects', { selector: 'h2' });
};

export const rowNames = () => {
  const names: string[] = [];
  for (const row of screen.queryAllByTestId('manage-projects-row')) {
    names.push(row.querySelector('[data-testid="manage-projects-name"]')?.textContent ?? '');
  }
  return names;
};

export const rowOf = (name: string) => {
  for (const row of screen.getAllByTestId('manage-projects-row')) {
    if (row.querySelector('[data-testid="manage-projects-name"]')?.textContent === name) {
      return row;
    }
  }
  throw new Error(`no row for ${name}`);
};

// The names in the delete confirm table, in order.
export const confirmNames = () => {
  const names: string[] = [];
  for (const row of screen.queryAllByTestId('delete-projects-row')) {
    names.push(row.querySelector('[data-testid="delete-projects-name"]')?.textContent ?? '');
  }
  return names;
};

// The delete confirm row of a project, as its whole text.
export const confirmRowText = (name: string) => {
  for (const row of screen.getAllByTestId('delete-projects-row')) {
    if (row.querySelector('[data-testid="delete-projects-name"]')?.textContent === name) {
      return row.textContent;
    }
  }
  throw new Error(`no confirm row for ${name}`);
};
