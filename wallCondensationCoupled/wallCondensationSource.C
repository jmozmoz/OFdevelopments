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

\*---------------------------------------------------------------------------*/

#include "wallCondensationSource.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"
#include "basicThermo.H"
#include "coupledPolyPatch.H"
#include "surfaceInterpolate.H"
#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(wallCondensationSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        wallCondensationSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::wallCondensationSource::wallCondensationSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(name, modelType, dict, mesh),
    specieName_(coeffs_.lookupOrDefault<word>("specie", "none")),
    filmMassSourceFluid_(
        IOobject
        (
            "filmMassSource" + specieName_,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "filmMassSource" + specieName_,
            dimDensity/dimTime,
            Zero
        )
    ),
    filmEnergySourceFluid_(
        IOobject
        (
            "filmEnergySource" + specieName_,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "filmEnergySource" + specieName_,
            dimEnergy/dimVolume/dimTime,
            Zero
        )
    )

{
    const basicThermo& thermo =
        mesh_.lookupObject<basicThermo>(basicThermo::dictName);

    fieldNames_.setSize(3);

    fieldNames_[0] = thermo.he().name();
    fieldNames_[1] = specieName_;
    fieldNames_[2] = "rho";
//    fieldNames_[3] = "p_rgh"; // we do not need to activate the source
//                              // for p_rgh, because it will be activated
//                              // automatically, if rho is active

    Info<<"Field names: " << fieldNames_ << endl;

    applied_.setSize(fieldNames_.size(), false);}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::wallCondensationSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (true)
    {
        Info<< type() << ": applying source 1 to " << eqn.psi().name()
            << " for fieldi: " << fieldi << endl;
    }

    eqn += filmMassSourceFluid_;
}

void Foam::fv::wallCondensationSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (true)
    {
        Info<< type() << ": applying source 2 to " << eqn.psi().name()
            << " for fieldi: " << fieldi << endl;
    }

    if (fieldi == 0)
    {
        eqn += filmEnergySourceFluid_;
    }
    else
    {
        eqn += filmMassSourceFluid_;
    }
}


bool Foam::fv::wallCondensationSource::read(const dictionary& dict)
{
    if (cellSetOption::read(dict))
    {
        coeffs_.lookup("specie") >> specieName_;
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
