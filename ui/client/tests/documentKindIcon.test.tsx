import { describe, it, expect, afterEach } from 'vitest';
import { render, cleanup } from '@testing-library/react';
import { DocumentKindIcon } from '../src/components/project/DocumentKindIcon';
import { TabKind } from '../src/shared/tabKind';

afterEach(() => { cleanup(); });

// The tree icons come from the same source as the tabs (TAB_KIND_STYLES), so this
// also guards that the tree and the tabs stay aligned.
const cases: [TabKind, string][] = [
  [TabKind.Agent, 'StorageIcon'],
  [TabKind.Workflow, 'AccountTreeIcon'],
  [TabKind.Notebook, 'MenuBookIcon'],
  [TabKind.ProjectConfig, 'SettingsIcon'],
  [TabKind.Document, 'DescriptionOutlinedIcon'],
  [TabKind.Repository, 'SourceIcon'],
  [TabKind.CentralRepository, 'FolderOpenIcon'],
  [TabKind.Toolkit, 'HandymanIcon'],
];

describe('DocumentKindIcon', () => {
  it.each(cases)('renders the shared tab icon for %s', (kind, testId) => {
    const { getByTestId } = render(<DocumentKindIcon kind={kind} />);
    expect(getByTestId(testId)).toBeTruthy();
  });
});
