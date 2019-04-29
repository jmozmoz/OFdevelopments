/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "Fick.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ThermoType>
void Foam::Fick<ThermoType>::updateCoefficients()
{
    DijModel_().update();

    forAll(species(), i)
    {
        if (i != inertIndex_)
        {
            volScalarField tmpGamma = 0 / Dij(0,0);
            forAll(species(), j)
            {
                if (j != i)
                {
                    tmpGamma += x_[j] / Dij(i,j);
                }
            }

            D_[i] = thermo_.rho() * (1-x_[i]) / (tmpGamma + dimensionedScalar("SMALL", dimensionSet(0,-2,1,0,0), Foam::SMALL));

            volScalarField& yi = thermo_.composition().Y(i);
            forAll(D_[i], cellI)
            {
                if((1-yi[cellI])==0)
                {
                    D_[i][cellI] = 0;
                }
            }
            forAll(D_[i].boundaryField(), boundaryI)
            {
                forAll(D_[i].boundaryField()[boundaryI], faceI)
                {
                    if((1-yi.boundaryField()[boundaryI][faceI])==0)
                    {
                        D_[i].boundaryFieldRef()[boundaryI][faceI] = 0;
                    }
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::Fick<ThermoType>::Fick
(
    // A.Alexiou 2014
    //hsCombustionThermo& thermo,
    psiReactionThermo& thermo,
    const compressible::turbulenceModel& turbulence,
    const surfaceScalarField& phi
)
:
    multiSpeciesTransportModel(thermo, turbulence, phi),

    speciesThermo_
    (
        dynamic_cast<const multiComponentMixture<ThermoType>&>
            (this->thermo_).speciesData()
    )

{
    D_.setSize(species().size()-1);

    updateMolarFractions();

    forAll(D_, i)
    {
        D_.set
        (
            i, new volScalarField
            (
                IOobject
                (
                    "D_" + species()[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("D", dimensionSet(1, -1, -1, 0, 0), 0.0)
            )
        );
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::scalar Foam::Fick<ThermoType>::correct
(
    fvMatrix<scalar>& yEqn,
    label fieldNumber
)
{
    scalar maxResidual = 0;
    scalar eqnResidual = 1;

    volScalarField& Yi = thermo_.composition().Y(fieldNumber); //yEqn.psi();

    yEqn -= fvm::laplacian(turbulence_.muEff(), Yi);
    yEqn += fvm::laplacian(D_[fieldNumber], Yi, "laplacian(D,Yi)");
    yEqn += Sy_[fieldNumber];

//        eqnResidual = yEqn.solve(mesh_.solver("Yi")).initialResidual();
//        maxResidual = max(eqnResidual, maxResidual);

//    ni = yEqn.flux();

//    nt -= ni;

//    n_[inertIndex_] = nt;

    return maxResidual;
}

template<class ThermoType>
bool Foam::Fick<ThermoType>::read()
{
    if (regIOobject::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
