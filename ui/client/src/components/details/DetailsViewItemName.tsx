import { RenameField } from '../../elements/RenameField';

export const DetailsViewItemName = ({
  itemKey,
  setItemKey = undefined,
}: {
  itemKey: string,
  setItemKey?: (newKey: string | undefined) => void | undefined,
}) => {
  return (
    <RenameField
      value={itemKey}
      setValue={setItemKey}
      labelMinWidth="100px"
    />
  )
}
