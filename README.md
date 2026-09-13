# Finite-Volume CFD and Heat Transfer Solvers

This repository collects numerical-method implementations developed from graduate coursework in computational fluid dynamics and heat transfer. The projects are organized by physical and numerical topic rather than by homework number.

Each completed project includes:

- the governing equation and boundary conditions;
- a concise description of the discretization or approximation;
- reproducible MATLAB or Python source code;
- a single command that regenerates the reported results;
- comparison against an analytical solution, correlation, or conservation check.

## Projects

### [1D Poisson Equations: Weighted Residuals and Finite Volumes](01-1d-poisson-weighted-residual-fvm/README.md)

Galerkin, point-collocation, and subdomain-collocation approximations are compared with analytical solutions for two source profiles. A finite-volume solution of a discontinuous-source problem is then verified against its piecewise analytical solution.

### [Finite-Volume Heat Conduction in a Composite Wall](02-composite-wall-conduction/README.md)

A conservative cell-centered formulation resolves volumetric heat generation, an insulated boundary, a material interface, and convection to the surroundings. Harmonic face conductance, TDMA, analytical verification, and global energy conservation are included.

### [Implicit Finite-Volume Simulation of Convective Slab Cooling](03-transient-slab-cooling/README.md)

A fully implicit transient solver predicts the cooling of a hot plane wall with convection at both surfaces. The implementation includes half control volumes at the boundaries, TDMA time marching, event interpolation, conservation checks, and comparison with the analytical eigenfunction solution.

### [SOR Solution of 2D Steady Heat Conduction with Mixed Boundary Conditions](04-2d-steady-conduction-sor/README.md)

A cell-centered finite-volume solver treats isothermal, insulated, and convective boundaries on a rectangular domain. SOR relaxation factors are benchmarked, grid refinement is compared with a separation-of-variables solution, and the global boundary heat balance is verified.

### [Transient Heating of a Concrete-Steel Composite with Convection and Radiation](05-composite-wall-fire-heating/README.md)

A fully implicit finite-volume solver models a concrete-steel composite exposed to a 900 degC environment. The material interface is resistance-based, radiation is handled through converged nonlinear iteration, and coupled grid/time refinement and energy conservation are verified.

### [Thermally Developing Laminar Channel Flow: Scheme and Solver Comparison](06-thermally-developing-channel-flow/README.md)

A two-dimensional finite-volume solution of the parallel-plate thermal-entry problem compares four convection schemes and three stationary iterative solvers, with grid refinement and a Shah-London local Nusselt-number benchmark.

### [Laminar Flat-Plate Flow and Heat Transfer: A Solver Revisited](07-flat-plate-flow-revisited/README.md)

A reconstructed staggered finite-volume SIMPLEC solver couples momentum, pressure and heat transport over a partially heated flat plate. The project documents errors in the original unsuccessful implementation and includes manufactured operator tests, conservation checks, grid and domain studies, and comparisons with Blasius and a heat-transfer correlation.

The original course file names, student identifiers, assignment PDFs, and third-party solution files are intentionally excluded from this public-facing structure.
