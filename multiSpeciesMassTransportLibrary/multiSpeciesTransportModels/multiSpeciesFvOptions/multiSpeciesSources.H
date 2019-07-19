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
    Foam::fv::multiSpeciesTransportSources

Description
    Corrections according to Kumar G et al. to account correctly
    calculate the energy evolution in case of multi-species mixtures
    with different Cp per specimen.

    Two terms are added to the energy equation:
    correct -fvm::laplacian(turbulence->alphaEff(), he)
     = -\nabla (\lapha_{eff} \nabla h)
    to only include diffusive heat conduction due to temperature
    gradient: -\nabla (\lambda_{eff} \nabla T)

    References:
    \verbatim
       V. Kumar G et al.
       Development and verification of a multi-species gas transport solver
       OpenFOAM Workshop 2019, Duisburg, Germany
    \endverbatim

Usage
    Example usage:
    \verbatim
    modifyDiffusion
    {
        type          multiSpeciesSources;
        selectionMode all;
    }
    \endverbatim

SourceFiles
    multiSpeciesSources.C

\*---------------------------------------------------------------------------*/

#ifndef multiSpeciesSources_H
#define multiSpeciesSources_H

#include "fvMesh.H"
#include "volFields.H"
#include "cellSetOption.H"
#include "basicThermo.H"
#include "basicSpecieMixture.H"
#include "multiSpeciesTransportModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                 Class multiSpeciesSources Declaration
\*---------------------------------------------------------------------------*/

class multiSpeciesSources
:
    public cellSetOption
{

private:

    // Private data

        const basicThermo& thermo_;

        volScalarField kappa_;

        DimensionedField<scalar, volMesh> chemistrySh_;

        autoPtr<multiSpeciesTransportModel> mstm_;

        const word thermoTypeThermo_;

        multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields_;


//
//        //- Current time index (used for updating)
//        label curTimeIndex_;

    // Private Member Functions

        //- Return the composition
        const basicSpecieMixture& composition() const;

        //- Update the model
        void update(const volScalarField& Cp);

        //- Helper function to apply to the energy equation
        template<class RhoFieldType>
        void apply(const RhoFieldType& rho, fvMatrix<scalar>& eqn);

        //- Disallow default bitwise copy construct
        multiSpeciesSources(const multiSpeciesSources&);

        //- Disallow default bitwise assignment
        void operator=(const multiSpeciesSources&);


public:

    //- Runtime type information
    TypeName("multiSpeciesTransportSources");


    // Constructors

        //- Construct from explicit source name and mesh
        multiSpeciesSources
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


            //- Correct the energy field
            virtual void correct(volScalarField& he);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
