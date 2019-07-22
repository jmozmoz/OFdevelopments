How to implement a new boundary condition and modify a solver without modifying
the solver source code
============

In this tutorial a step by step guide will be given how to implement a new
boundary condition for conjugate heat transfer simulations modelling the
condensation of vapor at a wall. In the implemented model the energy field as
well as the species fraction of vapor have to be modified. Vapor is removed
from the fluid and at the same time the released heat of condensation has to
be transferred to the fluid and the wall. Additionally, modifications of the
species fractions equation and the energy equation are necessary to correctly
simulate the diffusion of energy if the species have different heat capacities
to avoid unphysical temperatures.

Several aspects of the implementation are of general interest and can be used
for implementations of other boundary conditions. One is how to create a
boundary condition derived for the OpenFOAM class `mixedFvPatchScalarField`
used to model the energy balance at the surface between the fluid and the
solid domain. Another is how to modify the vapor fraction in the fluid domain
near the wall by source terms in different field equations. A way will be
shown how to implement this by using the `fvOptions` method, so no modifications
of the solver source code are necessary. An important point is to determine all
equations for which source terms have to be created. Another topic is how to
access all necessary fluid and wall fields and properties at the surface to
calculate the condensation rate and to store the source terms. This includes
the different ways to access values of mesh faces and cells like `patchField`,
`internalField`, `faceCells`.

The last topic discussed is how to use functions objects together with Python
and the libraries `pandas`, `numpy`, and `matplotlib` to postprocess the
results to verify and validate the model and its implementation.

The [slides](doc/tutorial-boundary-condition.pdf) are available in this
repository.