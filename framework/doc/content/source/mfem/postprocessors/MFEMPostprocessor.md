# MFEMPostprocessor

!if! function=hasCapability('mfem')

## Summary

Base class for MFEM postprocessors used to evaluate a single scalar.

## Overview

MFEM postprocessors calculate scalar quantities from the (aux)variables, typically after each
timestep.

An `MFEMPostprocessor` is derived from `MFEMExecutedObject`. Its ordering relative to MFEM initial
conditions, aux kernels, transfers, and other MFEM postprocessors is determined automatically from
detected data dependencies instead of manual execution groups.

The value of a postprocessor added to an MFEM problem is available as a spatially uniform scalar
coefficient named after it, and so may be used anywhere a scalar coefficient is expected. This
applies to postprocessors of all types, not only to those derived from `MFEMPostprocessor`, so that
values calculated elsewhere and transferred in, such as from a subapp into a
[Receiver.md] postprocessor, may also be used to build coefficients.

`MFEMPostprocessor` is a purely virtual base class. Derived classes
should override the `execute` and `getValue` methods.

!if-end!

!else
!include mfem/mfem_warning.md
