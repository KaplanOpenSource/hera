import { AccountTree, Code, Description, FolderOpen, Handyman, Science, Settings, Source } from '@mui/icons-material';
import { SvgIconComponent } from '@mui/icons-material';
import { TabKind } from './tabKind';

interface TabKindStyle {
  icon: SvgIconComponent;
  color: string;
  // Lighter variant used on the dark tab bar so the label/icon stays readable.
  darkColor: string;
  bg: string;
  bgOpacity?: number;
  selectedBgOpacity?: number;
}

export const TAB_KIND_STYLES: Record<TabKind, TabKindStyle> = {
  [TabKind.Notebook]: { icon: Code, color: '#4a6b3a', darkColor: '#9cc77f', bg: '126, 154, 110' },
  [TabKind.Document]: { icon: Description, color: '#3a5f80', darkColor: '#82b1d8', bg: '106, 140, 175' },
  [TabKind.Agent]: { icon: Science, color: '#7a4a76', darkColor: '#cf9fca', bg: '176, 122, 171' },
  [TabKind.Workflow]: { icon: AccountTree, color: '#2f6f73', darkColor: '#7fc9cd', bg: '94, 168, 176' },
  [TabKind.ProjectConfig]: { icon: Settings, color: '#555555', darkColor: '#b8b8b8', bg: '138, 138, 138', bgOpacity: 0.15, selectedBgOpacity: 0.25 },
  [TabKind.Repository]: { icon: Source, color: '#7a5530', darkColor: '#d0a877', bg: '192, 145, 94' },
  [TabKind.CentralRepository]: { icon: FolderOpen, color: '#6e5c30', darkColor: '#cdb87f', bg: '176, 151, 94' },
  [TabKind.Toolkit]: { icon: Handyman, color: '#4a7a6b', darkColor: '#8fcfbc', bg: '94, 176, 152' },
};

export const buildTabKindCss = (dark: boolean) => {
  return Object.entries(TAB_KIND_STYLES)
    .map(([kind, style]) => {
      const opacity = style.bgOpacity ?? 0.18;
      const selectedOpacity = style.selectedBgOpacity ?? 0.30;
      return `.tab-kind-${kind} {
  background-color: rgba(${style.bg}, ${opacity});
  color: ${dark ? style.darkColor : style.color};
}
.tab-kind-${kind}.flexlayout__tab_button--selected {
  background-color: rgba(${style.bg}, ${selectedOpacity});
}`;
    })
    .join('\n\n');
};
