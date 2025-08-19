# Table of Contents

1.  [What it is](#org8e27d27)
2.  [Installation](#org5380b22)
3.  [License](#orgf26070d)
4.  [Contributors](#org0388d7f)
5.  [Further reading](#org820f59e)
6.  [This document](#orgca16204)


<a id="org8e27d27"></a>

# What it is

The purpose of this library is to support working with the OpenSource
CFD-toolbox [OpenFOAM](http://www.openfoam.org) and its forks

This [Python](http://www.python.org) library can be used to

-   analyze the logs produced by OpenFoam-solvers
-   execute OpenFoam-solvers and utilities and analyze their output
    simultaneously
-   manipulate the parameter files and the initial-conditions of a run
    in a non-destructive manner
-   plots the residuals of OpenFOAM solvers
-   lots of other stuff

Most of this functionality is made available to the user in the form
of command-line utilities.

PyFoam does all this strictly "from the outside": by writing parameter
files and reading the output of the solvers. Without compiled parts or
being linked to OpenFOAM.

More information is found on [the OpenFOAM Wiki](http://openfoamwiki.net/index.php/Contrib_PyFoam).
Introductory presentations on PyFoam can be found there


<a id="org5380b22"></a>

# Installation

The easiest way to install PyFoam is the Python package-manager `pip`:

    pip install PyFoam

which will install PyFoam


<a id="orgf26070d"></a>

# License

PyFoam is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.  See the file COPYING in this directory,
for a description of the GNU General Public License terms under which
you can copy the files.


<a id="org0388d7f"></a>

# Contributors

If not otherwise noted in a source-files the primary author is Bernhard Gschaider

The people who contributed to PyFoam (If I forgot someone: tell me):

-   Bernhard Gschaider
-   Martin Beaudoin
-   Fabian Pollesböck
-   Etienne Lorriaux
-   Bruno Santos
-   Marc Immer
-   Oliver Borm


<a id="org820f59e"></a>

# Further reading

These documents give further information

-   **ReleaseNotes:** list of the changes between versions (newest
    versions are on top).
-   **DeveloperNotes:** document with information for people who want to
    contribute to `PyFoam`

For information on the usage see the presentations on [the `PyFOAM`
page on the OpenFOAM-Wiki](https://openfoamwiki.net/index.php/Contrib_PyFoam)


<a id="orgca16204"></a>

# This document

The original source of this document is `README.org`. The
`README.md` (to which the `README` links) is automatically generated
and should **not** be edited. The reason for this setup is that most
Web-GUIs for VCS insist on Markdown as a markup language


