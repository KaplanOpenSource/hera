import { KeyboardArrowDown, KeyboardArrowUp, KeyboardDoubleArrowDown } from '@mui/icons-material';
import { ReactNode } from 'react';
import { ButtonTooltip } from '../../elements/ButtonTooltip';

// How much of a document's detail section is shown above the editor:
//   Both — read-only header + the editable field tree
//   Tree — field tree only (header hidden)
//   None — both hidden, giving the editor the most space
export enum DetailsVisibility {
  Both = 'both',
  Tree = 'tree',
  None = 'none',
}

// Each press reveals less, then wraps back to showing everything.
const NEXT: { [key in DetailsVisibility]: DetailsVisibility } = {
  [DetailsVisibility.Both]: DetailsVisibility.Tree,
  [DetailsVisibility.Tree]: DetailsVisibility.None,
  [DetailsVisibility.None]: DetailsVisibility.Both,
};

// The chevron points down when more is visible, up when less; the double
// chevron marks the fully-expanded state.
const ICON: { [key in DetailsVisibility]: ReactNode } = {
  [DetailsVisibility.Both]: <KeyboardDoubleArrowDown />,
  [DetailsVisibility.Tree]: <KeyboardArrowDown />,
  [DetailsVisibility.None]: <KeyboardArrowUp />,
};

const TITLE: { [key in DetailsVisibility]: string } = {
  [DetailsVisibility.Both]: 'Hide header',
  [DetailsVisibility.Tree]: 'Hide fields',
  [DetailsVisibility.None]: 'Show details',
};

export const DetailsVisibilityToggle = ({
  value,
  onChange,
}: {
  value: DetailsVisibility,
  onChange: (value: DetailsVisibility) => void,
}) => {
  return (
    <ButtonTooltip
      title={TITLE[value]}
      onClick={() => onChange(NEXT[value])}
    >
      {ICON[value]}
    </ButtonTooltip>
  );
};
