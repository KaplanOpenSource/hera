import { Done } from '@mui/icons-material';
import { SimpleTreeView } from '@mui/x-tree-view';
import { ProjectDocument } from '@shared/types';
import { useState } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';
import { execPython } from '../../io/execPython';
import { DetailsViewItem } from './DetailsViewItem';

const FORBIDDEN_FIELDS = ['_id', '_cls', 'projectName'];

export const DetailsViewDocument = ({
  doc,
  setDoc,
}: {
  doc: any,
  setDoc: (newVal: any) => void,
}) => {
  const [shownDoc, setShownDoc] = useState<any>(JSON.parse(JSON.stringify(doc)));
  const name = doc?.desc?.datasourceName || doc?.type || doc._cls;

  const updateDocument = async () => {
    const docid = (doc as ProjectDocument)._id.$oid;
    const lines = [`
from hera.datalayer import All
doc = All.getDocumentByID('${docid}')
`];
    for (const [field, prevVal] of Object.entries(doc)) {
      if (!FORBIDDEN_FIELDS.includes(field) && JSON.stringify(prevVal) !== JSON.stringify(shownDoc[field])) {
        lines.push(`doc.${field} = ${JSON.stringify(shownDoc[field])}`)
      }
    }
    lines.push(`
doc.save()
docs = All.getDocumentByID('${docid}')
result = docs.asDict(with_id=True)
`)
    const code = lines.join('\n');
    const { data } = await execPython(code);
    if (data) {
      setDoc(data)
    }
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
