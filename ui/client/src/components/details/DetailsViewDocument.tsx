export const DetailsViewDocument = ({
  doc,
}: {
  doc: any,
}) => {
  return (<>
    <pre>{JSON.stringify(doc, null, 2)}</pre>
  </>)
}