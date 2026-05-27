import { Code, Description, FolderOpen, Handyman, Science, Settings, Source } from '@mui/icons-material';
import { SvgIconComponent } from '@mui/icons-material';
import { TabKind } from './tabKind';

interface TabKindStyle {
  icon: SvgIconComponent;
  color: string;
  bg: string;
  bgOpacity?: number;
  selectedBgOpacity?: number;
}

export const TAB_KIND_STYLES: Record<TabKind, TabKindStyle> = {
  [TabKind.Notebook]: { icon: Code, color: '#4a6b3a', bg: '126, 154, 110' },
  [TabKind.Document]: { icon: Description, color: '#3a5f80', bg: '106, 140, 175' },
  [TabKind.Agent]: { icon: Science, color: '#7a4a76', bg: '176, 122, 171' },
  [TabKind.ProjectConfig]: { icon: Settings, color: '#555555', bg: '138, 138, 138', bgOpacity: 0.15, selectedBgOpacity: 0.25 },
  [TabKind.Repository]: { icon: Source, color: '#7a5530', bg: '192, 145, 94' },
  [TabKind.CentralRepository]: { icon: FolderOpen, color: '#6e5c30', bg: '176, 151, 94' },
  [TabKind.Toolkit]: { icon: Handyman, color: '#4a7a6b', bg: '94, 176, 152' },
};

export const tabKindCss = Object.entries(TAB_KIND_STYLES)
  .map(([kind, style]) => {
    const opacity = style.bgOpacity ?? 0.18;
    const selectedOpacity = style.selectedBgOpacity ?? 0.30;
    return `.tab-kind-${kind} {
  background-color: rgba(${style.bg}, ${opacity});
  color: ${style.color};
}
.tab-kind-${kind}.flexlayout__tab_button--selected {
  background-color: rgba(${style.bg}, ${selectedOpacity});
}`;
  })
  .join('\n\n');
