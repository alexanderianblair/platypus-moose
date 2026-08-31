# MFEMIntegralConstraint

!if! function=hasCapability('mfem')

## Summary

Base class for objects weakly constraining a scalar integral quantity of an MFEM variable.

## Overview

Where [MFEMEssentialConstraint.md] imposes a constraint *strongly*, by eliminating degrees of
freedom from the discrete system, classes deriving from `MFEMIntegralConstraint` impose a single
scalar equation *weakly*, by introducing a multiplier. The multiplier is held by an
[MFEMScalarVariable.md] named by `scalar_variable`, and is solved for alongside the field
variables.

Writing the multiplier as $\lambda$ and the block of the system belonging to the constrained
variable as $A u = b$, the constrained system solved is

\begin{equation}
\begin{bmatrix} A & c \\ c^T & d \end{bmatrix}
\begin{bmatrix} u \\ \lambda \end{bmatrix} =
\begin{bmatrix} b \\ t \end{bmatrix}
\end{equation}

so that the multiplier adds $\lambda c$ to the residual of the constrained variable, and its own
row imposes

\begin{equation}
c^T u + d \lambda = t.
\end{equation}

Derived classes supply the coupling vector $c$ and the diagonal $d$; the target $t$ is the value
of the constrained quantity. Because the constraint row is the transpose of the coupling column,
the augmented system remains symmetric if the field blocks are.

A pure Lagrange multiplier constraint is the special case $d = 0$. A non-zero $d$ arises whenever
the multiplier is itself a physical unknown contributing to the constrained integral quantity, as
the loop voltage does for a closed conductor carrying a prescribed net current; see
[MFEMNetCurrentIntegralConstraint.md].

Essential degrees of freedom of the constrained variable are eliminated from the constraint in the
same way as from the rest of the system: their entries in $c$ are zeroed, and their known
contribution $c^T u$ is moved onto $t$.

## Limitations

The multiplier blocks are added to the assembled block matrix, so an integral constraint requires
`assembly_level = legacy` (the default). The augmented system is symmetric *indefinite*, so it
cannot be preconditioned with `MFEMHypreAMS` or `MFEMHypreBoomerAMG`; a direct solver such as
`MFEMSuperLU` should be used instead. Integral constraints are not currently supported in complex
(time-harmonic), eigenvalue, transient, or nonlinear MFEM problems, and each constraint requires
its own scalar variable.

Concrete implementations:

- [MFEMNetCurrentIntegralConstraint.md]

!if-end!

!else
!include mfem/mfem_warning.md
