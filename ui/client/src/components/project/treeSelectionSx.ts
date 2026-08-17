import { alpha, Theme } from '@mui/material';

// Selection styling shared by the documents and repositories trees: the selected
// row's name and kind icon turn cyan (the accent) over a subtle slate fill with a
// thin ring. The ring is an inset box-shadow so it doesn't shift the layout.
export const treeSelectionSx = (theme: Theme) => {
  const isDark = theme.palette.mode === 'dark';
  // In dark mode use the mockup's slate-navy selection; in light mode fall back to
  // a soft tint of the primary color.
  const selBg = isDark ? '#122336' : alpha(theme.palette.primary.main, 0.14);
  const selRing = isDark ? '#1e293b' : alpha(theme.palette.primary.main, 0.45);
  const focusRing = isDark ? '#334155' : alpha(theme.palette.primary.main, 0.6);
  return {
    '& .MuiTreeItem-content.Mui-selected': {
      backgroundColor: selBg,
      borderRadius: theme.shape.borderRadius,
      boxShadow: `inset 0 0 0 1px ${selRing}`,
      color: theme.palette.primary.main,
      '& .MuiSvgIcon-root': { color: theme.palette.primary.main },
    },
    '& .MuiTreeItem-content.Mui-selected:hover': {
      backgroundColor: selBg,
    },
    '& .MuiTreeItem-content.Mui-selected.Mui-focused': {
      backgroundColor: selBg,
      boxShadow: `inset 0 0 0 1px ${focusRing}`,
    },
    '& .MuiTreeItem-content.Mui-selected.Mui-focused:hover': {
      backgroundColor: selBg,
    },
  };
};
