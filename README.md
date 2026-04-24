# 1D Diffusion Equation Solver: Explicit & Implicit Finite Difference Methods

This repository contains a modular MATLAB script for numerically solving the 1D diffusion (heat conduction) equation. It implements both the **Explicit** and **Implicit** Finite Difference Methods (FDM) to explore transient thermodynamic processes, mass conservation, and numerical stability.

This code was developed as part of a Master's level Computational Materials Science laboratory to demonstrate the transition from conditionally stable explicit time-stepping to unconditionally stable implicit matrix inversion (via the Thomas Algorithm).

## 🚀 Features & The "Scenario Switchboard"

Instead of maintaining multiple scattered scripts, this project utilizes a **Scenario Switchboard** architecture. By changing a single variable (`scenario = X`) at the top of the script, the code automatically reconfigures the physics, boundary conditions, time steps, and visualization tools.

**Available Scenarios Include:**
* **Standard Evolution:** Simulates the diffusion of a Gaussian concentration profile over time.
* **Numerical Instability:** Deliberately violates the explicit Fourier stability limit ($\Delta t \le \Delta x^2 / 2D$) to demonstrate catastrophic exponential error amplification.
* **Implicit Unconditional Stability:** Proves the robustness of the implicit scheme by applying massive time steps that would crash an explicit solver.
* **Modified Initial Conditions (Block/Top-Hat):** Simulates the rapid smoothing of infinite concentration gradients (e.g., a localized precipitate in a pure matrix).
* **Modified Boundary Conditions:** Toggles between Zero-Flux (Neumann) and Fixed-Concentration (Dirichlet) boundaries to observe mass conservation versus steady-state gradient formation.

## 📊 Visualization
Depending on the selected scenario, the script automatically generates:
1. **Real-time 2D Animations:** Tracking the concentration profile dynamically across the spatial domain.
2. **3D Spatiotemporal Surface Plots:** Providing a holistic view of the entire processing history.
3. **Static 1x3 Subplots:** Formatted specifically for clean extraction into academic lab reports or publications.

## 🛠️ How to Use

1. Clone the repository and open the script in MATLAB (or MATLAB Online).
2. Locate the **SCENARIO SWITCHBOARD** at the very top of the script.
3. Uncomment exactly **one** scenario line (e.g., `scenario = 1;`).
4. Ensure all other scenario lines are commented out with `%`.
5. Run the script.

## ⚙️ Requirements
* MATLAB (Base installation). No additional toolboxes are strictly required.

## 🧠 Mathematical Core
* **Explicit Scheme:** Forward-Time Central-Space (FTCS) approximation.
* **Implicit Scheme:** Backward-Time Central-Space (BTCS) approximation, solved in $O(N)$ time complexity utilizing a custom implementation of the **Thomas Algorithm** for tridiagonal matrix systems.
