                     *multiSpeciesTransportModels README*

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
*                                                                             *
* SUMMARY                                                                     *
*                                                                             *
* 1) diffusivityModels                                                        *
*                                                                             *
* 2) multiSpeciesTransportModels                                              *
*                                                                             *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


-------------------------------------------------------------------------------

1) diffusivityModels

   - binaryDiffusivityModel -> base class for all binary
     diffusivity models

   - constant -> binary diffusivity model with constant
     diffusion coefficients

   - Fuller -> binary diffusivity model based on 
     Fuller-Schettler-Giddings correlation

   - ChapmanEnskog -> binary diffusivity model based on
     Champan-Enskog correlation

   - Wilke -> binary diffusivity model based on Wilke-Lee
     correlation

   - Knudsen -> Knudsen diffusivity model

   - diffusivityModel -> class that collects the binary
     diffusion coefficients for a set of species

   - KnudsenDiffusivityModel -> class that collects the
     Knudsen diffusion coefficients for a set of species


-------------------------------------------------------------------------------

2) multiSpeciesTransportModels

   - multiSpeciesTransportModels/multiSpeciesTransportModel -> base class for
     all multiSpecies transport models

   - multiSpeciesTransportModels/Fick -> diffusive mass fluxes are calculated
     using the common Fick law for multicomponent mixture (D_alpha is function
     of molar/mass fractions)

   - multiSpeciesTransportModels/FickDilutedMixture -> similar to Fick model,
     but D_alpha IS NOT function of molar/mass fractions

   - multiSpeciesTransportModels/SchmidtNumber -> diffusive mass fluxes are
     calculated using the Schmidt number

   - multiSpeciesTransportModels/Lewis -> diffusive mass fluxes are
     calculated using the Lewis number

   - multiSpeciesTransportModels/Bosanquet -> similar to Fick model,
     but D_alpha is corrected with Knudsen effects

   - multiSpeciesTransportModels/MaxwellStefan -> diffusive mass fluxes are
     calculated using the Maxwell-Stefan correlation (j_alpha depend by
     gradient of all species)
