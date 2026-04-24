# 1D Diffusion Equation Solver: Explicit & Implicit Finite Difference Methods

This repository contains a modular MATLAB script for numerically solving the 1D diffusion (heat conduction) equation. It implements both the **Explicit** and **Implicit** Finite Difference Methods (FDM) to explore transient thermodynamic processes, mass conservation, and numerical stability.

This code was developed as part of a Master's level Computational Materials Science laboratory to demonstrate the transition from conditionally stable explicit time-stepping to unconditionally stable implicit matrix inversion (via the Thomas Algorithm).

## 🧠 Mathematical Core
* **Explicit Scheme:** Forward-Time Central-Space (FTCS) approximation.
* **Implicit Scheme:** Backward-Time Central-Space (BTCS) approximation, solved in $O(N)$ time complexity utilizing a custom implementation of the **Thomas Algorithm** for tridiagonal matrix systems.
