import { Done } from '@mui/icons-material';
import { SimpleTreeView } from '@mui/x-tree-view';
import { ProjectDocument } from '@shared/types';
import { useState } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { DetailsViewItem } from './DetailsViewItem';
import { execPython } from '../../io/execPython';

export const DetailsViewDocument = ({
  doc,
}: {
  doc: any,
}) => {
  const [shownDoc, setShownDoc] = useState<any>(JSON.parse(JSON.stringify(doc)));
  const name = doc?.desc?.datasourceName || doc?.type || doc._cls;

  const updateDocument = async () => {
    const id = (doc as ProjectDocument)._id.$oid;
    const lines = [`
from hera.datalayer import All
doc = All.getDocumentByID('${id}')
`];
    const FORBIDDEN_FIELDS = ['_id', '_cls', 'projectName'];
    for (const [field, prevVal] of Object.entries(doc)) {
      if (!FORBIDDEN_FIELDS.includes(field) && JSON.stringify(prevVal) !== JSON.stringify(shownDoc[field])) {
        lines.push(`doc.${field} = ${JSON.stringify(shownDoc[field])}`)
      }
    }
    lines.push('doc.save()')
    const code = lines.join('\n');
    await execPython(code);
  }

  return (
    <SimpleTreeView>
      <DetailsViewItem
        itemKey={name}
        itemValue={shownDoc}
        level={0}
        index={0}
        setItemValue={newVal => setShownDoc(newVal)}
        labelAdditions={<>
          {JSON.stringify(doc) === JSON.stringify(shownDoc) ? null : (
            <ButtonTooltip
              title='Update Document'
              onClick={updateDocument}
            >
              <Done />
            </ButtonTooltip>
          )}
        </>}
      />
    </SimpleTreeView>
  )
}
