import { describe, it, expect, vi, beforeEach, afterEach } from 'vitest';
import { render, cleanup, waitFor } from '@testing-library/react';
import { MemoryRouter } from 'react-router-dom';

vi.mock('../src/io/fetchPython', () => ({
  fetchPython: vi.fn(),
  fetchPythonClean: vi.fn(),
}));

import { fetchPython } from '../src/io/fetchPython';
import { FetchProjects } from '../src/io/FetchProjects';
import { useProjectStore } from '../src/stores/useProjectStore';
import { useProjectListStore } from '../src/stores/useProjectListStore';

const renderFetch = () => {
  return render(
    <MemoryRouter initialEntries={['/Alpha']}>
      <FetchProjects urlProjectName="Alpha" />
    </MemoryRouter>
  );
};

beforeEach(() => {
  vi.mocked(fetchPython).mockReset();
  vi.mocked(fetchPython).mockResolvedValue({ data: null });
  useProjectListStore.setState({ projectNames: [{ name: 'Alpha' }, { name: 'Beta' }] });
});

afterEach(() => {
  cleanup();
});

describe('FetchProjects after startup', () => {
  it('does not fetch the project startup already loaded', async () => {
    useProjectStore.setState({
      currProjectName: 'Alpha',
      currProject: { name: 'Alpha', documents: [] },
    });

    renderFetch();

    await waitFor(() => {
      expect(useProjectStore.getState().currProjectName).toBe('Alpha');
    });
    expect(fetchPython).not.toHaveBeenCalled();
  });

  it('fetches when the store holds another project', async () => {
    useProjectStore.setState({
      currProjectName: 'Alpha',
      currProject: { name: 'Beta', documents: [] },
    });

    renderFetch();

    await waitFor(() => {
      expect(fetchPython).toHaveBeenCalled();
    });
    const labels = vi.mocked(fetchPython).mock.calls[0].map(c => c.label);
    expect(labels).toEqual(['toolkits', 'project Alpha']);
  });
});
