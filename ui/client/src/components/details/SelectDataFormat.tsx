import { SelectProperty } from "../../elements/SelectProperty";
import { useServerConstants } from "../../stores/useServerConstants";

export const SelectDataFormat = ({
  value,
  setValue,
}: {
  value: string,
  setValue: (v: string) => void,
}) => {
  const { dataTypes } = useServerConstants();
  return (
    <SelectProperty
      value={value}
      label="dataFormat"
      setValue={v => setValue(v)}
      menuItems={Object.entries(dataTypes).map(([_upcasename, name]) => ({ name }))}
    />
  )
}