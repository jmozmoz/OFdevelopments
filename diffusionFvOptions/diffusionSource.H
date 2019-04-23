/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::fv::diffusionSource

Description
    Corrections according to Kelm et al. to account correctly
    calculate the energy evolution in case of multi-species mixtures
    with different Cp per specimen.

    Two terms are added to the energy equation:
    correct -fvm::laplacian(turbulence->alphaEff(), he)
     = -\nabla (\lapha_{eff} \nabla h)
    to only include diffusive heat conduction due to temperature
    gradient: -\nabla (\lambda_{eff} \nabla T)

    References:
    \verbatim
        S. Kelm et al. (2019(
        Open Source CFD Modeling of Containment Transport and Mixing Processes
        Development and First Validation of ‘containmentFOAM’
        International CFD Network Meeting, GRS Garching, March 19th 2019
    \endverbatim

Usage
    Example usage:
    \verbatim
    modifyDiffusion
    {
        type       modifyDiffusion;
        selectionMode all;
    }
    \endverbatim

SourceFiles
    diffusionSource.C
    diffusionSourceIO.C

\*---------------------------------------------------------------------------*/

#ifndef diffusionSource_H
#define diffusionSource_H

#include "fvMesh.H"
#include "volFields.H"
#include "cellSetOption.H"
#include "NamedEnum.H"
#include "basicThermo.H"
#include "psiReactionThermo.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class diffusionSource Declaration
\*---------------------------------------------------------------------------*/

class diffusionSource
:
    public cellSetOption
{
public:

//    enum thermoMode
//    {
//        mdThermo,
//        mdLookup
//    };

//    static const NamedEnum<thermoMode, 2> thermoModeTypeNames_;


private:

    // Private data

        const basicThermo& thermo_;
        const psiReactionThermo& psiThermo_;
        const basicSpecieMixture& composition_;
        const PtrList<volScalarField>& Y_;
        tmp<volScalarField> gammaTmp_;
        const word thermoTypeThermo_;

//
//        //- Current time index (used for updating)
//        label curTimeIndex_;
//
//        //- Temperature change cached for source calculation when alpha1 updated
//        scalarField deltaT_;


    // Private Member Functions

        //- Update the model
        void update(const volScalarField& Cp);

        //- Helper function to apply to the energy equation
        template<class RhoFieldType>
        void apply(const RhoFieldType& rho, fvMatrix<scalar>& eqn);

        //- Disallow default bitwise copy construct
        diffusionSource(const diffusionSource&);

        //- Disallow default bitwise assignment
        void operator=(const diffusionSource&);


public:

    //- Runtime type information
    TypeName("modifyDiffusion");


    // Constructors

        //- Construct from explicit source name and mesh
        diffusionSource
        (
            const word& sourceName,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    // Member Functions

        // Add explicit and implicit contributions to compressible equation

            //- Add explicit contribution to compressible enthalpy and species
            // equation
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldi
            );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "diffusionSourceTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
