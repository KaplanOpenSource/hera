import { SvgIconProps } from '@mui/material';
import { TabKind } from '../../shared/tabKind';
import { TAB_KIND_STYLES } from '../../shared/tabKindConfig';

// Uses the same per-kind icon as the document tabs (TAB_KIND_STYLES), so the tree
// and the tabs always show matching icons.
export const DocumentKindIcon = ({
  kind,
  ...props
}: {
  kind: TabKind,
} & SvgIconProps) => {
  const Icon = TAB_KIND_STYLES[kind].icon;
  return <Icon {...props} />;
};
