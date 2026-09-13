# Development History and Authorship

## Origin of the work

I began these projects in ME613 Computational Fluid Dynamics. During the
course I worked through the governing equations, discretizations, boundary
conditions, MATLAB or Python implementations, parameter studies, plots, and
written reports. The original working directory contains the dated submission
reports and intermediate files, but those records are not redistributed in
this public repository.

The public repository is a later portfolio edition. I used my course codes and
reports to reconstruct the problem statements, then standardized filenames,
separated reusable solvers from plotting scripts, regenerated results, and
added explicit verification checks. This process also exposed mistakes that
were not recognized when the assignments were submitted.

## Project lineage

| Public project | Coursework origin | Portfolio revision |
|---|---|---|
| 01 - 1D Poisson methods | Weighted-residual and finite-volume assignments | Unified related exercises and checked them against analytical solutions |
| 02 - Composite-wall conduction | Steady one-dimensional conduction assignment | Corrected interface treatment and added convergence and conservation checks |
| 03 - Transient slab cooling | Implicit transient-conduction assignment | Added event interpolation, analytical comparison, and energy verification |
| 04 - 2D conduction with SOR | Mixed-boundary two-dimensional assignment | Added relaxation benchmarking, direct-solver comparison, and heat balance |
| 05 - Composite fire heating | Transient concrete-steel assignment | Reworked the material interface and nonlinear radiation treatment |
| 06 - Developing channel heat transfer | First term project | Corrected coefficient dimensions and node distances; reran scheme, solver, and grid comparisons |
| 07 - Heated flat-plate flow | Second term project, which did not converge | Diagnosed the original indexing and pressure-correction defects and built a new verified reconstruction |

## Later review and assistance

The portfolio cleanup and the Project 07 reconstruction were completed with
AI-assisted code review and implementation support. I retained responsibility
for selecting the coursework to present, providing the original source and
reports, deciding the physical scope, reviewing the numerical interpretation,
and accepting the final public version. The verification scripts are included
so that numerical claims can be reproduced rather than accepted from the
README alone.

Project 07 is deliberately described as a later reconstruction. Its converged
results are not presented as results obtained in the original course
submission. A third-party solution notebook was consulted to resolve an
ambiguous heating-start location and reference correlation; it is credited in
that project's README and is not included as project code.
