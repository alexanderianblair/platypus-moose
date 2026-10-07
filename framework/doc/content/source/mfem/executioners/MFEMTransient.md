# MFEMTransient

!if! function=hasCapability('mfem')

## Overview

`MFEMTransient` is the `Executioner` class used to solve time dependent MFEM finite element
problems, calling the [MFEMProblemSolve.md] solve object to execute one or more
MFEM `TimeDependentProblemOperators`.

As in all `Executioner` classes using the [MFEMProblemSolve.md] solve object,
the desired device and assembly level to use during problem set-up and solution can be selected.

The time integration scheme is selected with either the `scheme` parameter, shared with other
transient executioners, or the `mfem_scheme` parameter, which provides further schemes implemented by
MFEM. The schemes available through `scheme` map onto MFEM's ODE solvers as follows:

| `scheme` | MFEM ODE solver |
| - | - |
| `implicit-euler` (default) | `BackwardEulerSolver` |
| `explicit-euler` | `ForwardEulerSolver` |
| `crank-nicolson` | `TrapezoidalRuleSolver` |
| `explicit-midpoint` | `RK2Solver` with $a = 1/2$ |
| `dirk` | `SDIRK23Solver`, L-stable and second order |
| `explicit-tvd-rk-2` | `RK2Solver` with $a = 1$ (Heun's method) |

The schemes available through `mfem_scheme` are the implicit midpoint rule, singly diagonally
implicit Runge-Kutta (SDIRK and ESDIRK) schemes, the generalized-alpha scheme, whose damping of
high frequencies is set by `rho_inf`, explicit Runge-Kutta schemes of orders three and four, and
the implicit-explicit (IMEX) forward-backward Euler scheme. The IMEX scheme treats kernels and
integrated boundary conditions explicitly if their `implicit` parameter is `false`, and implicitly
otherwise. The `implicit` parameter is ignored by all other schemes, which treat every term either
implicitly or explicitly.

Schemes with explicit stages, including the first stage of `crank-nicolson`, invert the mass
operator, so require every equation to contain a time derivative of its trial variable. See
[TimeDependentEquationSystem.md] for the form of the stages solved.

The [TimeStepper](/TimeSteppers/index.md) system is fully supported providing wide choice for the timestep(s) `dt`.

## Example Input File Syntax

!listing test/tests/mfem/kernels/heattransfer.i block=Executioner

The convergence of each scheme is verified on a problem with time-dependent essential and
convective boundary conditions:

!listing test/tests/mfem/time_integrators/mms.i

!syntax parameters /Executioner/MFEMTransient

!syntax inputs /Executioner/MFEMTransient

!syntax children /Executioner/MFEMTransient

!if-end!

!else
!include mfem/mfem_warning.md
