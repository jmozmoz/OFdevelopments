/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2018 OpenFOAM Foundation
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
    Foam::fv::wallCondensationSource

Description


Usage
    Example usage:
    \verbatim
    wallCondensationSource1
    {
        type            wallCondensationSource;
        active          yes;
        selectionMode   all;
        specie          H2O;
    }
    \endverbatim

Note
    structure taken from effectivenessHeatExchangerSource
    /OpenFOAM-6/src/fvOptions/sources/derived/effectivenessHeatExchangerSource/effectivenessHeatExchangerSource.H

SourceFiles
    wallCondensationSource.C

\*---------------------------------------------------------------------------*/

#ifndef wallCondensationSource_H
#define wallCondensationSource_H

#include "cellSetOption.H"
#include "autoPtr.H"
#include "interpolation2DTable.H"
#include "volFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
              Class wallCondensationSource Declaration
\*---------------------------------------------------------------------------*/

class wallCondensationSource
:
    public cellSetOption
{

protected:

    // Protected data

        //- Name of the affected specie
        word specieName_;

        //- Mass source term for specie removed from fluid
        volScalarField filmMassSourceFluid_;

        //- Energy source term for specific energy of specie removed from fluid
        volScalarField filmEnergySourceFluid_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        wallCondensationSource
        (
            const wallCondensationSource&
        );

        //- Disallow default bitwise assignment
        void operator=(const wallCondensationSource&);


public:

    //- Runtime type information
    TypeName("wallCondensationSource");


    // Constructors

        //- Construct from components
        wallCondensationSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~wallCondensationSource()
    {}


    // Member Functions

        // Explicit and implicit source

            //- Scalar
            virtual void addSup
            (
                fvMatrix<scalar>& eqn,
                const label fieldi
            )
            {
                NotImplemented;
            }


        // Explicit and implicit source for compressible equation

            //- Scalar
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldi
            );


        // IO

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
