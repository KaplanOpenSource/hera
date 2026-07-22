import { styled, Switch } from '@mui/material';
import { ThemeMode } from '../stores/useViewSettingsStore';

// Material Design filled sun / moon glyphs, drawn white inside the thumb.
const sunPath =
  'M12 7c-2.76 0-5 2.24-5 5s2.24 5 5 5 5-2.24 5-5-2.24-5-5-5zM2 13h2c.55 0 1-.45 1-1s-.45-1-1-1H2c-.55 0-1 .45-1 1s.45 1 1 1zm18 0h2c.55 0 1-.45 1-1s-.45-1-1-1h-2c-.55 0-1 .45-1 1s.45 1 1 1zM11 2v2c0 .55.45 1 1 1s1-.45 1-1V2c0-.55-.45-1-1-1s-1 .45-1 1zm0 18v2c0 .55.45 1 1 1s1-.45 1-1v-2c0-.55-.45-1-1-1s-1 .45-1 1zM5.99 4.58c-.39-.39-1.03-.39-1.41 0-.39.39-.39 1.03 0 1.41l1.06 1.06c.39.39 1.03.39 1.41 0s.39-1.03 0-1.41L5.99 4.58zm12.37 12.37c-.39-.39-1.03-.39-1.41 0-.39.39-.39 1.03 0 1.41l1.06 1.06c.39.39 1.03.39 1.41 0 .39-.39.39-1.03 0-1.41l-1.06-1.06zm1.06-10.96c.39-.39.39-1.03 0-1.41-.39-.39-1.03-.39-1.41 0l-1.06 1.06c-.39.39-.39 1.03 0 1.41s1.03.39 1.41 0l1.06-1.06zM7.05 18.36c.39-.39.39-1.03 0-1.41-.39-.39-1.03-.39-1.41 0l-1.06 1.06c-.39.39-.39 1.03 0 1.41s1.03.39 1.41 0l1.06-1.06z';
const moonPath =
  'M12 3c-4.97 0-9 4.03-9 9s4.03 9 9 9 9-4.03 9-9c0-.46-.04-.92-.1-1.36-.98 1.37-2.58 2.26-4.4 2.26-2.98 0-5.4-2.42-5.4-5.4 0-1.81.89-3.42 2.26-4.4-.44-.06-.9-.1-1.36-.1z';

const glyph = (path: string) =>
  `url('data:image/svg+xml;utf8,<svg xmlns="http://www.w3.org/2000/svg" height="18" width="18" viewBox="0 0 24 24"><path fill="%23fff" d="${path}"/></svg>')`;

// A switch whose thumb is a sun in light mode and a moon in dark mode.
const SunMoonSwitch = styled(Switch)(() => ({
  width: 62,
  height: 34,
  padding: 7,
  '& .MuiSwitch-switchBase': {
    margin: 1,
    padding: 0,
    transform: 'translateX(6px)',
    '&.Mui-checked': {
      transform: 'translateX(28px)',
      '& .MuiSwitch-thumb': {
        backgroundColor: '#cbd5e1',
        '&:before': { backgroundImage: glyph(moonPath) },
      },
      '& + .MuiSwitch-track': {
        opacity: 1,
        backgroundColor: '#3a4a63',
      },
    },
  },
  '& .MuiSwitch-thumb': {
    backgroundColor: '#f9c22e',
    width: 32,
    height: 32,
    '&:before': {
      content: "''",
      position: 'absolute',
      width: '100%',
      height: '100%',
      left: 0,
      top: 0,
      backgroundRepeat: 'no-repeat',
      backgroundPosition: 'center',
      backgroundImage: glyph(sunPath),
    },
  },
  '& .MuiSwitch-track': {
    opacity: 1,
    backgroundColor: '#aab4c2',
    borderRadius: 20 / 2,
  },
}));

export const ThemeModeSwitch = ({
  mode,
  setMode,
}: {
  mode: ThemeMode,
  setMode: (mode: ThemeMode) => void,
}) => {
  return (
    <SunMoonSwitch
      checked={mode === ThemeMode.Dark}
      onChange={(e) => {
        e.stopPropagation();
        setMode(e.target.checked ? ThemeMode.Dark : ThemeMode.Light);
      }}
      slotProps={{ input: { 'aria-label': 'toggle dark mode' } }}
    />
  );
};
