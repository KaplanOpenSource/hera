import {
  AccountTree,
  DescriptionOutlined,
  InsertDriveFileOutlined,
  MenuBook,
  Settings,
  Storage,
} from '@mui/icons-material';
import { SvgIconProps } from '@mui/material';
import { ElementType } from 'react';
import { TabKind } from '../../shared/tabKind';

// A distinct MUI icon per document kind shown in the tree.
const KIND_ICON: { [key in TabKind]?: ElementType } = {
  [TabKind.Agent]: Storage,
  [TabKind.Workflow]: AccountTree,
  [TabKind.Notebook]: MenuBook,
  [TabKind.ProjectConfig]: Settings,
  [TabKind.Document]: DescriptionOutlined,
};

export const DocumentKindIcon = ({
  kind,
  ...props
}: {
  kind: TabKind,
} & SvgIconProps) => {
  const Icon = KIND_ICON[kind] ?? InsertDriveFileOutlined;
  return <Icon {...props} />;
};
