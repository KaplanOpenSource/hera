import { describe, it, expect, vi, afterEach } from 'vitest';
import { render, screen, fireEvent, cleanup } from '@testing-library/react';
import { RawViewToggle } from '../src/components/details/RawViewToggle';

afterEach(() => { cleanup(); });

const switchInput = (container: HTMLElement) =>
  container.querySelector('input[type="checkbox"]') as HTMLInputElement;

describe('RawViewToggle', () => {
  it('renders a "Raw View" labeled switch, unchecked by default', () => {
    const { container } = render(<RawViewToggle rawView={false} setRawView={() => {}} />);
    expect(screen.getByText('Raw View')).toBeTruthy();
    expect(switchInput(container).checked).toBe(false);
  });

  it('shows checked when rawView is true', () => {
    const { container } = render(<RawViewToggle rawView={true} setRawView={() => {}} />);
    expect(switchInput(container).checked).toBe(true);
  });

  it('calls setRawView with the new value when toggled on', () => {
    const setRawView = vi.fn();
    const { container } = render(<RawViewToggle rawView={false} setRawView={setRawView} />);
    fireEvent.click(switchInput(container));
    expect(setRawView).toHaveBeenCalledWith(true);
  });

  it('calls setRawView(false) when toggled off', () => {
    const setRawView = vi.fn();
    const { container } = render(<RawViewToggle rawView={true} setRawView={setRawView} />);
    fireEvent.click(switchInput(container));
    expect(setRawView).toHaveBeenCalledWith(false);
  });
});
