// Python that scans the project's files directory for .ipynb files and
// registers any that are not already tracked as notebook documents.
export const buildDetectNotebooksCode = ({
  projectName,
  filesDir,
}: {
  projectName: string,
  filesDir: string,
}) => {
  return `
import os, glob
from hera.datalayer import All, Cache_Collection
filesDir = os.path.expanduser('${filesDir}')
existing = set()
for doc in Cache_Collection().getDocuments(projectName='${projectName}', type='notebook'):
    if doc.resource:
        existing.add(os.path.abspath(os.path.expanduser(doc.resource)))
found = glob.glob(os.path.join(filesDir, '**', '*.ipynb'), recursive=True)
found = [f for f in found if '.ipynb_checkpoints' not in f]
addedCount = 0
for f in found:
    if os.path.abspath(f) in existing:
        continue
    name = os.path.splitext(os.path.basename(f))[0]
    Cache_Collection().addDocument(
        projectName='${projectName}',
        resource=f,
        dataFormat="JSON_dict",
        type="notebook",
        desc={"datasourceName": name},
    )
    addedCount += 1
docs = All.getDocumentsAsDict('${projectName}', with_id=True)
project = {"name": '${projectName}', "documents": docs['documents']}
`;
};
