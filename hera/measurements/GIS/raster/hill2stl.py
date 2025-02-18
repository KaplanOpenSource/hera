#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 10:14:22 2025

@author: nirb

We create a STL file from a mathematical function.
We might want to change the code to allow it to use other parts of the repository e.g. create STL in topography
"""
import numpy as np

def function(x, y):
    """Define your function here."""
    return np.cos(x)**2 * np.cos(y)**2 + 3

def generate_solid_stl(f, x_range=(-20, 20), y_range=(-12, 12), resolution=100, filename="solid_output.stl"):
    """Generates a solid STL file with a flat bottom."""
    x = np.linspace(x_range[0], x_range[1], resolution)
    y = np.linspace(y_range[0], y_range[1], resolution)
    X, Y = np.meshgrid(x, y)
    Z = f(X, Y)  # Compute function values

    z_min = np.min(Z) - 0.5  # Set a flat base slightly below the surface

    # Open file for writing ASCII STL
    with open(filename, "w") as stl_file:
        stl_file.write("solid solid_surface\n")

        # Create top surface
        for i in range(resolution - 1):
            for j in range(resolution - 1):
                v1 = (X[i, j], Y[i, j], Z[i, j])
                v2 = (X[i, j+1], Y[i, j+1], Z[i, j+1])
                v3 = (X[i+1, j], Y[i+1, j], Z[i+1, j])
                v4 = (X[i+1, j+1], Y[i+1, j+1], Z[i+1, j+1])

                write_triangle(stl_file, v1, v2, v3)
                write_triangle(stl_file, v2, v4, v3)

        # Create flat base
        for i in range(resolution - 1):
            for j in range(resolution - 1):
                v1 = (X[i, j], Y[i, j], z_min)
                v2 = (X[i, j+1], Y[i, j+1], z_min)
                v3 = (X[i+1, j], Y[i+1, j], z_min)
                v4 = (X[i+1, j+1], Y[i+1, j+1], z_min)

                write_triangle(stl_file, v1, v3, v2)
                write_triangle(stl_file, v2, v3, v4)

        # Create vertical side walls
        for i in range(resolution - 1):
            for j in [0, resolution - 1]:  # Front and back
                v1 = (X[i, j], Y[i, j], Z[i, j])
                v2 = (X[i+1, j], Y[i+1, j], Z[i+1, j])
                v3 = (X[i, j], Y[i, j], z_min)
                v4 = (X[i+1, j], Y[i+1, j], z_min)

                write_triangle(stl_file, v1, v2, v3)
                write_triangle(stl_file, v2, v4, v3)

        for j in range(resolution - 1):
            for i in [0, resolution - 1]:  # Left and right
                v1 = (X[i, j], Y[i, j], Z[i, j])
                v2 = (X[i, j+1], Y[i, j+1], Z[i, j+1])
                v3 = (X[i, j], Y[i, j], z_min)
                v4 = (X[i, j+1], Y[i, j+1], z_min)

                write_triangle(stl_file, v1, v2, v3)
                write_triangle(stl_file, v2, v4, v3)

        stl_file.write("endsolid solid_surface\n")
    
    print(f"STL file saved as {filename}")

def write_triangle(file, v1, v2, v3):
    """Writes a triangle to the STL file in ASCII format."""
    normal = compute_normal(v1, v2, v3)
    file.write(f"  facet normal {normal[0]:.6f} {normal[1]:.6f} {normal[2]:.6f}\n")
    file.write("    outer loop\n")
    file.write(f"      vertex {v1[0]:.6f} {v1[1]:.6f} {v1[2]:.6f}\n")
    file.write(f"      vertex {v2[0]:.6f} {v2[1]:.6f} {v2[2]:.6f}\n")
    file.write(f"      vertex {v3[0]:.6f} {v3[1]:.6f} {v3[2]:.6f}\n")
    file.write("    endloop\n")
    file.write("  endfacet\n")

def compute_normal(v1, v2, v3):
    """Computes the normal of a triangle given three vertices."""
    v1, v2, v3 = np.array(v1), np.array(v2), np.array(v3)
    normal = np.cross(v2 - v1, v3 - v1)
    normal = normal / np.linalg.norm(normal)  # Normalize
    return normal

# Run the function
minx=-2
maxx=2
miny=-3
maxy=3
filename='test1.stl'
generate_solid_stl(function, x_range=(minx, maxx), y_range=(miny, maxy), resolution=100, filename=filename)
