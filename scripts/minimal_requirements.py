#!/usr/bin/env python3
"""Print the requirement lines for the packages a chosen list names.

The chosen list holds names only, so a version is pinned in one place.
A name the full file does not have comes out as the chosen list wrote it.
Asking for the full file prints every pin in it.

    minimal_requirements.py requirements.txt requirements-minimal.txt
"""
import re
import sys


def package_name(line):
    return re.sub(r"[-_.]+", "-", re.split(r"[\[=<>!;~ ]", line.strip())[0]).lower()


def minimal_requirements(full, chosen):
    pins = open(full).read().splitlines()
    for wanted in open(chosen):
        name = package_name(wanted)
        if name and name[0] not in "#-":
            found = [pin for pin in pins if package_name(pin) == name]
            if found:
                print(found[0])
            else:
                print(wanted.strip())


if __name__ == "__main__":
    minimal_requirements(*sys.argv[1:])
