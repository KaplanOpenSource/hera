import { render, screen, fireEvent, within, act } from '@testing-library/react';
import { ManageProjectsButton } from '../../src/components/header/ManageProjectsButton';

// Shared by the manage projects tests: the folders answer, opening the dialog, and finding rows.

export const FOLDERS = {
  projects: [
    { name: 'defaultProject', filesDirectory: '/data/default', folderExists: true },
    { name: 'Alpha', filesDirectory: '/data/alpha', folderExists: true },
    { name: 'Beta', filesDirectory: '/data/beta', folderExists: true },
    { name: 'Gamma', filesDirectory: '/data/gamma', folderExists: true },
  ],
};

export const openDialog = async () => {
  render(<ManageProjectsButton />);
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
