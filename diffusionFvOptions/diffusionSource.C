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

\*---------------------------------------------------------------------------*/

#include "diffusionSource.H"
#include "fvMatrices.H"
#include "fvmLaplacian.H"
#include "fvcLaplacian.H"
#include "fvcGrad.H"
#include "fvcDiv.H"
#include "fvmDdt.H"
#include "addToRunTimeSelectionTable.H"
#include "turbulentFluidThermoModel.H"
#include "psiReactionThermo.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{

    namespace fv
    {
        defineTypeNameAndDebug(diffusionSource, 0);

        addToRunTimeSelectionTable
        (
            option,
            diffusionSource,
            dictionary
        );
    }
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// copied from OpenFOAM-6/src/semiPermeableBaffle/derivedFvPatchFields/semiPermeableBaffleVelocity/semiPermeableBaffleVelocityFvPatchVectorField.C
const Foam::basicSpecieMixture&
Foam::fv::diffusionSource::composition() const
{
    const word& name = basicThermo::dictName;

    if (mesh_.foundObject<psiReactionThermo>(name))
    {
        return mesh_.lookupObject<psiReactionThermo>(name).composition();
    }
    else if (mesh_.foundObject<rhoReactionThermo>(name))
    {
        return mesh_.lookupObject<rhoReactionThermo>(name).composition();
    }
    else
    {
        FatalErrorInFunction
            << "Could not find a multi-component thermodynamic model."
            << exit(FatalError);

        return NullObjectRef<basicSpecieMixture>();
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::diffusionSource::diffusionSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    cellSetOption(sourceName, modelType, dict, mesh),
    thermo_(mesh_.lookupObject<basicThermo>(basicThermo::dictName)),
    gammaTmp_
    (
        new volScalarField(
            IOobject("gammaTmp", mesh_.time().timeName(), mesh_),
            mesh_,
            dimensionedScalar
            (
                "gammaTmp", dimMass/dimLength/dimTime, 0)
        )
    ),
    thermoTypeThermo_
    (
        mesh_.lookupObject<basicThermo>
        (
            basicThermo::dictName
        ).subDict
        (
            "thermoType"
        ).lookup("thermo")
    )
//    curTimeIndex_(-1),
//    deltaT_(cells_.size(), 0)
{

    const word inertSpecie(thermo_.lookup("inertSpecie"));
    const label inertIndex = composition().species()[inertSpecie];

    const PtrList<volScalarField>& Y = composition().Y();
    fieldNames_.setSize(Y.size());

    label j = 0;
    forAll(Y, i)
    {
        if (i != inertIndex) {
            fieldNames_[j++] = Y[i].name();
        }
    }

    fieldNames_[Y.size() - 1] = thermo_.he().name();

    Info<<"Field names: " << fieldNames_ << endl;

    applied_.setSize(fieldNames_.size(), false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::fv::diffusionSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (debug)
    {
        Info<< type() << ": applying source to " << eqn.psi().name()
            << " for fieldi: " << fieldi << endl;
    }

    const compressible::turbulenceModel& turbulence =
        mesh_.lookupObject<compressible::turbulenceModel>
        (
            turbulenceModel::propertiesName
        );

    if (fieldi == 0)
    {
        volScalarField& rho = const_cast<volScalarField&>(
            mesh_.lookupObject<volScalarField>("rho")
        );

        if (true)
        {
            const surfaceScalarField& phi =
                mesh_.lookupObject<surfaceScalarField>("phi");


            fvScalarMatrix rhoEqn
            (
                fvm::ddt(rho)
              + fvc::div(phi)
//              ==
//                fvOptions(rho)
            );

//            fvOptions.constrain(rhoEqn);

            rhoEqn.solve();

//            fvOptions.correct(rho);
        }
    }
    else if (0 < fieldi && fieldi < fieldNames_.size() - 1)
    {
//        const volScalarField& Yi = eqn.psi();
//        eqn -= fvm::laplacian(turbulence.muEff(), Yi);
    }
    else
    {
        const volScalarField& he = eqn.psi();
        if (thermoTypeThermo_ == "hConst"){
            // if Cp is temperature independent, the product rule can be applied
            // to \nabla he and the correction is>
            eqn -= fvc::laplacian
            (
                turbulence.alphaEff()*thermo_.T(), thermo_.Cp()
            );
        }
        else
        {
            // if Cp is temperature dpeendent, it is necessary to replace
            // laplacian alphaEff * he with laplacian lambdaEff T in the
            // energy equation
            eqn += fvc::laplacian
            (
                turbulence.alphaEff()*thermo_.Cp(), thermo_.T()
            );
            eqn -= fvc::laplacian(turbulence.alphaEff(), he);
        }

        const PtrList<volScalarField>& Y = composition().Y();

        forAll (Y, i)
        {
            tmp<volScalarField> hsTmp
            (
                new volScalarField(
                    IOobject("hsTmp", mesh_.time().timeName(), mesh_),
                    mesh_,
                    dimensionedScalar("hsTmp", dimEnergy/dimMass, 0)
                )
            );
            volScalarField& hs = hsTmp.ref();
            const volScalarField& p = mesh_.lookupObject<volScalarField>("p");
            forAll(Y[i], cellI)
            {
                hs[cellI] = composition().Hs(i, p[cellI], thermo_.T()[cellI]);
            }

            // the face values must be set, too. Otherwise diffusion
            // will be calculated correctly, if fixed values are set
            // on boundaries
            volScalarField::Boundary& hsBf = hs.boundaryFieldRef();

            forAll(hsBf, patchi)
            {
                scalarField& hsPatch = hsBf[patchi];
                const scalarField& pp = p.boundaryField()[patchi];
                const scalarField& Tp = thermo_.T().boundaryField()[patchi];

                forAll(hsPatch, facei)
                {
                    hsPatch[facei] = composition().Hs(i, pp[facei], Tp[facei]);
                }
            }


            eqn += fvc::laplacian(turbulence.muEff()*hs, Y[i]);
        }
    }
}


// ************************************************************************* //
