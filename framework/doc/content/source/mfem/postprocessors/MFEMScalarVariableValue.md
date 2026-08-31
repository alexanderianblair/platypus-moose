# MFEMScalarVariableValue

!if! function=hasCapability('mfem')

## Summary

Returns the value of an [MFEMScalarVariable.md].

## Overview

`MFEMScalarVariableValue` reports the single degree of freedom held by the named
[MFEMScalarVariable.md], so that a global unknown such as the multiplier of an
[MFEMIntegralConstraint.md] can be output or consumed by other MOOSE objects. The value is the
same on every rank, so no reduction is performed.

## Example Input File Syntax

!listing test/tests/mfem/constraints/closed_coil_net_current.i block=Postprocessors

!syntax parameters /Postprocessors/MFEMScalarVariableValue

!syntax inputs /Postprocessors/MFEMScalarVariableValue

!syntax children /Postprocessors/MFEMScalarVariableValue

!if-end!

!else
!include mfem/mfem_warning.md
