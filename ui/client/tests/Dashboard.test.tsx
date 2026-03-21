import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, screen, cleanup } from '@testing-library/react';
import { MemoryRouter } from 'react-router-dom';
import { Dashboard } from '../src/Dashboard';
import { useProjectStore } from '../src/stores/useProjectStore';

vi.mock('../src/io/FetchProjects', () => ({
  FetchProjects: () => null,
}));
vi.mock('../src/stores/useServerConstants', () => ({
  ServerConstantReader: () => null,
}));
vi.mock('../src/components/header/PageTitle', () => ({
  PageTitle: () => <div>PageTitle</div>,
}));
vi.mock('../src/components/header/ProjectChooser', () => ({
  ProjectChooser: () => <div>ProjectChooser</div>,
}));
vi.mock('../src/components/header/CommitIdShower', () => ({
  CommitIdShower: () => null,
}));
vi.mock('../src/components/project/ProjectTreeView', () => ({
  ProjectTreeView: () => <div>ProjectTreeView</div>,
}));
vi.mock('../src/components/details/DetailsViewPanel', () => ({
  DetailsViewPanel: () => <div>DetailsViewPanel</div>,
}));

const renderDashboard = () => {
  return render(
    <MemoryRouter>
      <Dashboard />
    </MemoryRouter>
  );
};

describe('Dashboard', () => {
  afterEach(() => {
    cleanup();
  });

  beforeEach(() => {
    useProjectStore.setState({
      projectNames: [],
      currProjectName: '* NONE *',
      currProject: null,
      toolkits: [],
    });
  });

  it('shows "No project loaded" when no project is selected', () => {
    renderDashboard();
    expect(screen.getByText('No project loaded')).toBeDefined();
  });

  it('shows project views when a project is loaded', () => {
    useProjectStore.setState({
      currProject: { name: 'TestProject', documents: [] },
    });
    renderDashboard();
    expect(screen.getByText('ProjectTreeView')).toBeDefined();
    expect(screen.getByText('DetailsViewPanel')).toBeDefined();
    expect(screen.queryByText('No project loaded')).toBeNull();
  });

  it('renders header components regardless of project state', () => {
    renderDashboard();
    expect(screen.getByText('PageTitle')).toBeDefined();
    expect(screen.getByText('ProjectChooser')).toBeDefined();
  });

  it('renders header components when a project is loaded', () => {
    useProjectStore.setState({
      currProject: { name: 'TestProject', documents: [] },
    });
    renderDashboard();
    expect(screen.getByText('PageTitle')).toBeDefined();
    expect(screen.getByText('ProjectChooser')).toBeDefined();
  });
});
