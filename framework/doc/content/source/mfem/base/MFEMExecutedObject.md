# MFEMExecutedObject

!if! function=hasCapability('mfem')

## Overview

`MFEMExecutedObject` is the base class for MFEM objects that are scheduled and executed by
[MFEMProblem.md] during a solve. It derives from [MFEMObject.md], [SetupInterface.md], and
[DependencyResolverInterface.md].

## Execution ordering

Execution order among `MFEMExecutedObject` instances is determined automatically by
[DependencyResolverInterface.md]: the MFEM scheduler builds a dependency graph from the variable,
postprocessor, and vector postprocessor names that each object declares it consumes or produces,
then performs a topological sort.

Derived class authors register dependency-bearing parameters at `validParams` time using the
provided helpers:

- `addDependencyParam<T>` — adds an optional parameter whose value names a consumed resource.
- `addRequiredDependencyParam<T>` — adds a required parameter whose value names a consumed resource.

The template parameter `T` must be one of `VariableName`, `std::vector<VariableName>`,
`PostprocessorName`, `std::vector<PostprocessorName>`, `VectorPostprocessorName`, or
`std::vector<VectorPostprocessorName>`. The scheduler infers the resource category from the type.

MFEM variables are also available as coefficients under their own names. Any scalar or vector
coefficient parameter whose value names an MFEM variable is therefore treated as consuming that
variable, without needing to be registered.

The resource supplied by this object is declared by overriding one of:

- `suppliedVariableName()` — name of the MFEM variable written by this object.
- `suppliedPostprocessorName()` — name of the postprocessor written by this object.
- `suppliedVectorPostprocessorName()` — name of the vector postprocessor written by this object.

Each returns `std::optional<std::string>`; return `std::nullopt` (the default) if the object supplies no resource in that category.

## Lifecycle methods

Each `MFEMExecutedObject` has three lifecycle methods that `MFEMProblem` calls in order for every
execution pass:

| Method | Purpose |
| - | - |
| `initialize()` | Pre-execution setup (zero accumulators, clear state, etc.). |
| `execute()` | Main work: read variables/coefficients, write results. |
| `finalize()` | Post-execution cleanup or reduction (e.g., parallel reduction for postprocessors). |

All three have empty default implementations. Derived classes override whichever they need.

## Derived classes

The following MFEM object types derive from `MFEMExecutedObject`:

- [MFEMAuxKernel.md] — updates real auxiliary variables.
- [MFEMComplexAuxKernel.md] — updates complex auxiliary variables.
- [MFEMPostprocessor.md] — computes a scalar quantity.
- [MFEMVectorPostprocessor.md] — computes an array of values.
- [MFEMSubMeshTransfer.md] — transfers variable data between a mesh and a submesh.
- `MFEMInitialCondition` — sets the initial value on an MFEM variable.

!if-end!

!else
!include mfem/mfem_warning.md
