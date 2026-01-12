export const idDocId = (oid: string) => {
  return `document_${oid}`;
}

export const idFromDocId = (docid: string) => {
  return docid && docid.startsWith("document_") ? docid.replace("document_", "") : undefined;
}
