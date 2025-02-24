# PolyDiM

**PolyDiM** (POLYtopal DIscretization Methods ) is a numerical computational library designed for solving partial differential equations (PDEs) using discretization methods that operate on generic polytopal (polygonal/polyhedral) meshes.

PolyDiM is inspired by and built upon the foundational research of the Numerical Analysis Group in the Department of Mathematical Sciences "Giuseppe Luigi Lagrange" (DISMA) at the Politecnico di Torino.

## Methods

The following methods have been implemented and tested:

- **FEM** (Finite Element Method) - A widely used numerical method for solving partial differential equations (PDEs), which discretizes a continuous domain into smaller local subdomains (finite elements).
- **VEM** (Virtual Element Method) - An advanced numerical method that extends FEM to general polygonal and polyhedral meshes, offering greater flexibility and high accuracy.

## Features

The library's key features include:

- **High-Order Approximation**: supports arbitrary polynomial degrees for highly accurate solutions.
- **General Mesh Support**: efficiently handles both polygonal and polyhedral meshes.
- **Modular and Extensible**: a well-structured framework enables users to define custom problems, mesh types, and solvers.
- **Applications**: suitable for solving linear and nonlinear elliptic PDEs, including Poisson, Elasticity, and Navier-Stokes problems.