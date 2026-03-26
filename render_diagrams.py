#!/usr/bin/env python3
"""
Pre-render all Mermaid diagrams in docs/ to SVG images using Docker.

Usage:
    python render_diagrams.py          # render all diagrams
    python render_diagrams.py --check  # just list diagrams without rendering

This script:
1. Finds all ```mermaid blocks in docs/**/*.md files
2. Renders each to SVG locally using the mermaid-cli Docker image (no internet needed)
3. Saves SVGs to docs/images/diagrams/
4. Replaces the mermaid code blocks with <img> tags pointing to the SVGs
5. Keeps the original mermaid source as an HTML comment for future edits

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
MERMAID_PATTERN = re.compile(r'```mermaid\n(.*?)```', re.DOTALL)
DOCKER_IMAGE = "minlag/mermaid-cli"


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

            input_file.write_text(mermaid_source)

            result = subprocess.run(
                [
                    "docker", "run", "--rm",
                    "-v", f"{tmpdir}:/data",
                    DOCKER_IMAGE,
                    "-i", "/data/input.mmd",
                    "-o", "/data/output.svg",
                    "-b", "transparent",
                    "-w", "2048",
                ],
                capture_output=True, text=True, timeout=30
            )

            if result.returncode != 0:
                print(f"  ERROR: {result.stderr.strip()[:200]}")
                return False

            if output_file.exists():
                # Copy to final destination
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


def process_file(md_path: Path, dry_run: bool = False) -> int:
    """Process a markdown file, replacing mermaid blocks with images."""
    with open(md_path, "r") as f:
        content = f.read()

    matches = list(MERMAID_PATTERN.finditer(content))
    if not matches:
        return 0

    rel_path = md_path.relative_to(DOCS_DIR)
    file_stem = str(rel_path).replace("/", "_").replace(".md", "")
    print(f"\n{rel_path} ({len(matches)} diagrams)")

    if dry_run:
        return len(matches)

    md_dir = md_path.parent
    rel_to_images = os.path.relpath(IMAGES_DIR, md_dir)

    replacements = []
    for i, match in enumerate(matches):
        mermaid_source = match.group(1).strip()
        diagram_hash = hashlib.md5(mermaid_source.encode()).hexdigest()[:8]
        img_name = f"{file_stem}_{i}_{diagram_hash}"
        img_path = IMAGES_DIR / f"{img_name}.svg"

        if render_mermaid_to_svg(mermaid_source, img_path):
            img_rel = f"{rel_to_images}/{img_name}.svg"
            replacement = (
                f'![Diagram]({img_rel})\n\n'
                f'<!-- mermaid source (for editing, paste into mermaid.live):\n'
                f'```mermaid\n{mermaid_source}\n```\n-->'
            )
            replacements.append((match.start(), match.end(), replacement))

    # Apply replacements in reverse order to preserve offsets
    for start, end, replacement in reversed(replacements):
        content = content[:start] + replacement + content[end:]

    with open(md_path, "w") as f:
        f.write(content)

    return len(replacements)


def main():
    dry_run = "--check" in sys.argv

    if not dry_run:
        check_docker()
        IMAGES_DIR.mkdir(parents=True, exist_ok=True)

    total = 0
    for md_path in sorted(DOCS_DIR.rglob("*.md")):
        count = process_file(md_path, dry_run=dry_run)
        total += count

    print(f"\n{'='*50}")
    if dry_run:
        print(f"Found {total} diagrams across docs/")
        print(f"Run without --check to render them.")
    else:
        print(f"Total: {total} diagrams rendered to {IMAGES_DIR}")
        print(f"\nTo edit a diagram:")
        print(f"  1. Find the <!-- mermaid source --> comment in the .md file")
        print(f"  2. Copy the mermaid code to https://mermaid.live")
        print(f"  3. Edit and re-run: python render_diagrams.py")


if __name__ == "__main__":
    main()
