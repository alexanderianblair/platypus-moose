# AddMFEMTimeIntegratorAction

!if! function=hasCapability('mfem')

## Overview

Action called to add an MFEM time integrator object to an MFEM problem, parsing content inside a
[`TimeIntegrators`](syntax/TimeIntegrators/index.md) block in the user input. The object it adds
selects the Runge-Kutta scheme used by the [MFEMTransient.md] executioner to advance the problem.
Only has an effect if the `Problem` type is set to [MFEMProblem.md].

## Example Input File Syntax

!listing test/tests/mfem/timeintegrators/mms.i block=Problem TimeIntegrators Executioner

!syntax parameters /TimeIntegrators/AddMFEMTimeIntegratorAction

!if-end!

!else
!include mfem/mfem_warning.md
