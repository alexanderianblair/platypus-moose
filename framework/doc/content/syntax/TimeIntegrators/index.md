# TimeIntegrators

!if! function=hasCapability('mfem')

The top-level `TimeIntegrators` syntax selects the time integration scheme used to advance a
transient MFEM problem. It is distinct from the `Executioner/TimeIntegrators` syntax, which
applies to libMesh-backed problems.

A transient MFEM problem is advanced by a single Runge-Kutta scheme, so the block takes exactly
one child block, whose name is the object name and whose `type` selects the scheme:

```text
[TimeIntegrators]
  [ti]
    type = MFEMRungeKuttaTimeIntegrator
    scheme = sdirk33
  []
[]
```

An arbitrary scheme may instead be given by its Butcher tableau with
[MFEMButcherTableauTimeIntegrator.md]:

```text
[TimeIntegrators]
  [ti]
    type = MFEMButcherTableauTimeIntegrator
    a = '0.25 0.0
         0.5  0.25'
    b = '0.5 0.5'
  []
[]
```

If no `TimeIntegrators` block is given, the problem is advanced with backwards Euler. Because the
scheme is then selected by the executioner's `scheme` parameter, setting that parameter and
providing a `TimeIntegrators` block at the same time is an error.

Only explicit and diagonally implicit tableaus are supported; a fully implicit scheme would
couple all stages into a single solve. See [MFEMTransient.md] for how a step is taken and
[TimeDependentEquationSystem.md] for the system solved at each stage.

!syntax list /TimeIntegrators objects=True actions=False subsystems=False

!if-end!

!else
!include mfem/mfem_warning.md
