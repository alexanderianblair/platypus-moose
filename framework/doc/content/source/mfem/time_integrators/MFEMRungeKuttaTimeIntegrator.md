# MFEMRungeKuttaTimeIntegrator

!if! function=hasCapability('mfem')

## Overview

Advances a transient MFEM problem with one of a set of named explicit or diagonally implicit
Runge-Kutta schemes, selected with [!param](/TimeIntegrators/MFEMRungeKuttaTimeIntegrator/scheme).
A scheme that is not listed here may be supplied by its Butcher tableau instead, with
[MFEMButcherTableauTimeIntegrator.md].

| `scheme` | Order | Stages | Properties |
| - | - | - | - |
| `implicit_euler` | 1 | 1 | L-stable, stiffly accurate |
| `implicit_midpoint` | 2 | 1 | A-stable, not stiffly accurate |
| `crank_nicolson` | 2 | 2 | A-stable, stiffly accurate, explicit first stage |
| `sdirk22` | 2 | 2 | L-stable, stiffly accurate |
| `sdirk33` | 3 | 3 | L-stable, stiffly accurate |
| `esdirk32` | 2 | 3 | L-stable, stiffly accurate, explicit first stage |
| `explicit_euler` | 1 | 1 | Explicit, conditionally stable |
| `rk4` | 4 | 4 | Explicit, conditionally stable |

Schemes that are stiffly accurate take their final state from the last stage, so essential
boundary values are satisfied exactly at the end of the timestep. Schemes that are not, such as
`implicit_midpoint`, combine the stage slopes with the quadrature weights instead, and their
essential boundary values at the end of the timestep are therefore extrapolated from the stage
times rather than imposed.

Every scheme listed here has stage order one, so on a stiff problem the observed order of
convergence is reduced towards two as the timestep grows relative to the fastest mode of the
spatial discretisation. The third order `sdirk33` in particular converges at its design order
only once the timestep resolves that mode.

Fully explicit schemes, and the explicit first stage of `crank_nicolson` and `esdirk32`, evaluate
their stage slope by solving a mass system. That requires a linear solver, and requires every
equation of the system to carry a time derivative kernel, since an equation without one has a
singular mass operator. The two schemes with only an explicit *first* stage remain L-stable, but
`explicit_euler` and `rk4` are conditionally stable, so their timestep is restricted by the
stiffest mode of the discretisation.

!syntax parameters /TimeIntegrators/MFEMRungeKuttaTimeIntegrator

!syntax inputs /TimeIntegrators/MFEMRungeKuttaTimeIntegrator

!syntax children /TimeIntegrators/MFEMRungeKuttaTimeIntegrator

!if-end!

!else
!include mfem/mfem_warning.md
