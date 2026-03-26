#!/usr/bin/env python3
"""
Pre-render Mermaid diagrams in docs/ to SVG images using Docker.

Usage:
    python render_diagrams.py          # render new/changed diagrams only
    python render_diagrams.py --force  # re-render all diagrams
    python render_diagrams.py --check  # list diagrams without rendering

Handles two formats:
  1. Fresh ```mermaid blocks (first run)
  2. Previously rendered blocks stored as HTML comments (re-render on change)

The content hash is embedded in the SVG filename, so unchanged diagrams
are skipped automatically (~0 seconds). Only changed diagrams trigger
Docker (~2 seconds each).

Requirements:
    Docker (the mermaid-cli image is pulled automatically)
    Run: make mermaid-pull   (or docker pull minlag/mermaid-cli)
"""

import os
import re
import sys
import hashlib
import subprocess
import tempfile
from pathlib import Path

DOCS_DIR = Path(__file__).parent / "docs"
IMAGES_DIR = DOCS_DIR / "images" / "diagrams"
DOCKER_IMAGE = "minlag/mermaid-cli"

# Pattern 1: fresh ```mermaid blocks (not yet rendered)
FRESH_PATTERN = re.compile(r'```mermaid\n(.*?)```', re.DOTALL)

# Pattern 2: previously rendered (image + HTML comment with mermaid source)
RENDERED_PATTERN = re.compile(
    r'!\[Diagram\]\((.*?\.svg)\)\n\n'
    r'<!-- mermaid source \(for editing, paste into mermaid\.live\):\n'
    r'```mermaid\n(.*?)```\n-->',
    re.DOTALL
)


def check_docker():
    """Verify Docker is available and the mermaid-cli image exists."""
    try:
        subprocess.run(["docker", "info"], capture_output=True, check=True)
    except (FileNotFoundError, subprocess.CalledProcessError):
        print("ERROR: Docker is not running. Start Docker and try again.")
        sys.exit(1)

    result = subprocess.run(
        ["docker", "images", "-q", DOCKER_IMAGE],
        capture_output=True, text=True
    )
    if not result.stdout.strip():
        print(f"Pulling {DOCKER_IMAGE}...")
        subprocess.run(["docker", "pull", DOCKER_IMAGE], check=True)


def render_mermaid_to_svg(mermaid_source: str, output_path: Path) -> bool:
    """Render a Mermaid diagram to SVG using Docker locally."""
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            input_file = Path(tmpdir) / "input.mmd"
            output_file = Path(tmpdir) / "output.svg"

            os.chmod(tmpdir, 0o777)
            input_file.write_text(mermaid_source)
            os.chmod(input_file, 0o666)

            uid = os.getuid()
            gid = os.getgid()

            result = subprocess.run(
                [
                    "docker", "run", "--rm",
                    "--user", f"{uid}:{gid}",
                    "-v", f"{tmpdir}:/data",
                    DOCKER_IMAGE,
                    "-i", "/data/input.mmd",
                    "-o", "/data/output.svg",
                    "-b", "transparent",
                    "-w", "2048",
                ],
                capture_output=True, text=True, timeout=60
            )

            if result.returncode != 0:
                print(f"  ERROR: {result.stderr.strip()[:200]}")
                return False

            if output_file.exists():
                output_path.parent.mkdir(parents=True, exist_ok=True)
                output_path.write_bytes(output_file.read_bytes())
                print(f"  Rendered: {output_path.name}")
                return True
            else:
                print(f"  ERROR: output file not created")
                return False

    except subprocess.TimeoutExpired:
        print(f"  ERROR: render timed out")
        return False
    except Exception as e:
        print(f"  ERROR: {e}")
        return False


def make_img_name(file_stem: str, index: int, mermaid_source: str) -> str:
    """Generate a deterministic image name from file + index + content hash."""
    diagram_hash = hashlib.md5(mermaid_source.encode()).hexdigest()[:8]
    return f"{file_stem}_{index}_{diagram_hash}"


def build_replacement(img_rel: str, mermaid_source: str) -> str:
    """Build the image + HTML comment replacement string."""
    return (
        f'![Diagram]({img_rel})\n\n'
        f'<!-- mermaid source (for editing, paste into mermaid.live):\n'
        f'```mermaid\n{mermaid_source}\n```\n-->'
    )


def process_file(md_path: Path, dry_run: bool = False, force: bool = False) -> tuple:
    """Process a markdown file. Returns (rendered_count, skipped_count)."""
    with open(md_path, "r") as f:
        content = f.read()

    rel_path = md_path.relative_to(DOCS_DIR)
    file_stem = str(rel_path).replace("/", "_").replace(".md", "")
    md_dir = md_path.parent
    rel_to_images = os.path.relpath(IMAGES_DIR, md_dir)

    # Collect all diagrams: previously rendered first, then fresh
    diagrams = []  # list of (start, end, mermaid_source, current_img_name_or_None)
    rendered_ranges = set()

    # First pass: find previously rendered blocks (image + HTML comment)
    for match in RENDERED_PATTERN.finditer(content):
        old_img_ref = match.group(1)
        mermaid_source = match.group(2).strip()
        old_img_name = Path(old_img_ref).stem
        diagrams.append((match.start(), match.end(), mermaid_source, old_img_name))
        rendered_ranges.add((match.start(), match.end()))

    # Second pass: find fresh ```mermaid blocks NOT inside a rendered block
    for match in FRESH_PATTERN.finditer(content):
        # Skip if this match falls inside an already-rendered block
        inside_rendered = any(
            rs <= match.start() and match.end() <= re
            for rs, re in rendered_ranges
        )
        if not inside_rendered:
            diagrams.append((match.start(), match.end(), match.group(1).strip(), None))

    # Sort by position
    diagrams.sort(key=lambda x: x[0])

    if not diagrams:
        return (0, 0)

    print(f"\n{rel_path} ({len(diagrams)} diagrams)")

    if dry_run:
        for i, (_, _, src, old_name) in enumerate(diagrams):
            new_name = make_img_name(file_stem, i, src)
            status = "unchanged" if old_name == new_name else "CHANGED" if old_name else "NEW"
            print(f"  [{status}] {new_name}")
        return (len(diagrams), 0)

    rendered = 0
    skipped = 0
    replacements = []

    for i, (start, end, mermaid_source, old_img_name) in enumerate(diagrams):
        new_img_name = make_img_name(file_stem, i, mermaid_source)
        img_path = IMAGES_DIR / f"{new_img_name}.svg"
        img_rel = f"{rel_to_images}/{new_img_name}.svg"

        if img_path.exists() and not force and old_img_name == new_img_name:
            print(f"  Skipped (unchanged): {new_img_name}.svg")
            skipped += 1
            # Still need to ensure the markdown has the right reference
            replacements.append((start, end, build_replacement(img_rel, mermaid_source)))
            continue

        # Clean up old SVG if hash changed
        if old_img_name and old_img_name != new_img_name:
            old_path = IMAGES_DIR / f"{old_img_name}.svg"
            if old_path.exists():
                old_path.unlink()
                print(f"  Removed old: {old_img_name}.svg")

        if render_mermaid_to_svg(mermaid_source, img_path):
            rendered += 1
            replacements.append((start, end, build_replacement(img_rel, mermaid_source)))

    # Apply replacements in reverse order
    for start, end, replacement in reversed(replacements):
        content = content[:start] + replacement + content[end:]

    with open(md_path, "w") as f:
        f.write(content)

    return (rendered, skipped)


def main():
    dry_run = "--check" in sys.argv
    force = "--force" in sys.argv

    if not dry_run:
        check_docker()
        IMAGES_DIR.mkdir(parents=True, exist_ok=True)

    total_rendered = 0
    total_skipped = 0
    for md_path in sorted(DOCS_DIR.rglob("*.md")):
        rendered, skipped = process_file(md_path, dry_run=dry_run, force=force)
        total_rendered += rendered
        total_skipped += skipped

    print(f"\n{'='*50}")
    if dry_run:
        print(f"Found {total_rendered} diagrams across docs/")
        print(f"Run without --check to render them.")
    else:
        print(f"Rendered: {total_rendered}, Skipped (unchanged): {total_skipped}")
        if total_rendered == 0 and total_skipped > 0:
            print(f"All diagrams are up to date.")
        print(f"\nTo edit a diagram:")
        print(f"  1. Edit the mermaid source in the <!-- comment --> in the .md file")
        print(f"  2. Re-run: python render_diagrams.py")
        print(f"  3. Only changed diagrams will be re-rendered (~2s each)")


if __name__ == "__main__":
    main()
