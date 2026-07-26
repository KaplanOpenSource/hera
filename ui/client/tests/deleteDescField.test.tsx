import { describe, it, expect, vi, afterEach } from 'vitest';
import { render, screen, fireEvent, cleanup } from '@testing-library/react';
import { DetailsViewDocumentContent } from '../src/components/details/DetailsViewDocumentContent';
import { DetailsViewDocument } from '../src/components/details/DetailsViewDocument';
import { ProjectObj, DocumentObj } from '../src/objects/ProjectObj';
import { copyOnly } from '../src/utils/utils';

const makeDoc = (desc: any) => {
  const data = { _id: { $oid: '1' }, _cls: 'C', type: 'sometype', desc, resource: '/r' } as any;
  const project = new ProjectObj({ name: 'P', documents: [data] } as any);
  return new DocumentObj(data, project);
};

const clickDeleteFor = (label: string) => {
  // The delete icon lives in the same tree row as the field name.
  const nameNode = screen.getByText(label);
  const row = nameNode.closest('li') ?? nameNode.parentElement!;
  const icon = row.querySelector('[data-testid="DeleteIcon"]');
  fireEvent.click(icon!.closest('button')!);
};

afterEach(() => cleanup());

describe('copyOnly', () => {
  it('keeps only the listed fields that exist', () => {
    expect(copyOnly({ a: 1, b: 2, c: 3 }, ['a', 'c', 'z'])).toEqual({ a: 1, c: 3 });
  });
});

describe('deleting a field under desc', () => {
  it('removes the field from desc while keeping hidden desc fields', () => {
    // toolkit is a HIDE_ON_DESC field: hidden from the tree in formulated view.
    const doc = makeDoc({ myField: 'hello', toolkit: 'GIS' });
    const setShownDoc = vi.fn();
    render(
      <DetailsViewDocumentContent
        doc={doc}
        setDoc={vi.fn()}
        shownDoc={doc.data}
        setShownDoc={setShownDoc}
      />
    );

    clickDeleteFor('myField');

    expect(setShownDoc).toHaveBeenCalledTimes(1);
    const newDoc = setShownDoc.mock.calls[0][0];
    // myField deleted, hidden toolkit preserved.
    expect(newDoc.desc).toEqual({ toolkit: 'GIS' });
  });

  it('does not offer a type chip on desc (its hidden fields make a switch unsafe)', () => {
    const doc = makeDoc({ myField: 'hello', toolkit: 'GIS' });
    render(
      <DetailsViewDocumentContent
        doc={doc}
        setDoc={vi.fn()}
        shownDoc={doc.data}
        setShownDoc={vi.fn()}
      />
    );
    // desc is the only top-level object; with no chip on it, "object" never shows.
    expect(screen.queryByText('object')).toBeNull();
    // but scalar fields still have a type chip
    expect(screen.getAllByText('string').length).toBeGreaterThan(0);
  });

  it('makes the row disappear in the live view', () => {
    const doc = makeDoc({ myField: 'hello', toolkit: 'GIS' });
    render(<DetailsViewDocument doc={doc} setDoc={vi.fn()} />);

    expect(screen.getByText('myField')).toBeDefined();
    clickDeleteFor('myField');
    expect(screen.queryByText('myField')).toBeNull();
  });
});
