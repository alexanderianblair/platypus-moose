# MFEMButcherTableauTimeIntegrator

!if! function=hasCapability('mfem')

## Overview

Advances a transient MFEM problem with the Runge-Kutta scheme given by the Butcher tableau
supplied in the input file. For an $s$ stage scheme,

!equation
k_i = f\left(t + c_i \delta t,\, u_n + \delta t \sum_j a_{ij} k_j\right), \qquad
u_{n+1} = u_n + \delta t \sum_i b_i k_i

where [!param](/TimeIntegrators/MFEMButcherTableauTimeIntegrator/a) holds the row-major entries
of the $s \times s$ matrix $a$, [!param](/TimeIntegrators/MFEMButcherTableauTimeIntegrator/b) the
quadrature weights, and [!param](/TimeIntegrators/MFEMButcherTableauTimeIntegrator/c) the stage
times as fractions of the timestep. If `c` is omitted it is taken to be the row sums of `a`.

Only explicit and diagonally implicit tableaus are supported, so every entry of `a` above the
diagonal must be zero; a fully implicit scheme would couple all $s$ stages into a single solve.
The first order consistency conditions, that `b` sums to one and that each entry of `c` matches
the corresponding row sum of `a`, are checked when the tableau is constructed.

A stage with a nonzero diagonal coefficient $a_{ii}$ is solved implicitly for the stage state.
A stage with a vanishing diagonal coefficient is evaluated by solving a mass system for the stage
slope, which requires every equation of the system to carry a time derivative kernel.

For example, the two stage, second order L-stable scheme of Alexander, which is also available as
`scheme = sdirk22` in [MFEMRungeKuttaTimeIntegrator.md], is

```text
[TimeIntegrators]
  [ti]
    type = MFEMButcherTableauTimeIntegrator
    a = '0.29289321881345248 0.0
         0.70710678118654752 0.29289321881345248'
    b = '0.70710678118654752 0.29289321881345248'
    c = '0.29289321881345248 1.0'
  []
[]
```

!syntax parameters /TimeIntegrators/MFEMButcherTableauTimeIntegrator

!syntax inputs /TimeIntegrators/MFEMButcherTableauTimeIntegrator

!syntax children /TimeIntegrators/MFEMButcherTableauTimeIntegrator

!if-end!

!else
!include mfem/mfem_warning.md
