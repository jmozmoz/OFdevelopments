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

#include "MaxwellStefan.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class ThermoType>
void Foam::MaxwellStefan<ThermoType>::updateCoefficients()
{
    DijModel_().update();

    const label N = species().size()-1;

    forAll(thermo_.composition().Y(0), cellI)
    {
        scalarRectangularMatrix A(N,N);
        scalarRectangularMatrix B(N,N);
        scalarRectangularMatrix C(N,N);

        for(label i=0; i<N; i++)
        {
            for(label j=0; j<N; j++)
            {
                if(j == i)
                {
                    A[i][j] =
                        -x_[i][cellI]/Dij(i,N)[cellI]/W(N);

                    for(label l=0; l<N+1; l++)
                    {
                        if(l != i)
                        {
                            A[i][j] -=
                                x_[l][cellI]/Dij(i,l)[cellI]/W(i);
                        }
                    }

                    B[i][j] =
                        -(x_[i][cellI]/W(N) + (1-x_[i][cellI])/W(i));
                }
                else
                {
                    A[i][j] = x_[i][cellI]
                        * (1/Dij(i,j)[cellI]/W(j) - 1/Dij(i,N)[cellI]/W(N));

                    B[i][j] = x_[i][cellI] * (1/W(j) - 1/W(N));
                }
            }
        }

        scalarRectangularMatrix invA = SVDinv(A);
        multiply(C,invA, B);
        scalarRectangularMatrix invC = SVDinv(C);

        // Prevent spurious numerical oscillations due to the matrix inversion
        for(label i=0; i<N; i++)
        {
            for(label j=0; j<N; j++)
            {
                D(i,j)[cellI] = thermo_.rho()()[cellI]*C[i][j];
                G(i,j)[cellI] = invC[i][j]/thermo_.rho()()[cellI];
            }
        }
    }

    forAll(thermo_.composition().Y(0).boundaryField(), patchI)
    {
        forAll(thermo_.composition().Y(0).boundaryField()[patchI], faceI)
        {
            scalarRectangularMatrix A(N,N);
            scalarRectangularMatrix B(N,N);
            scalarRectangularMatrix C(N,N);

            for(label i=0; i<N; i++)
            {
                for(label j=0; j<N; j++)
                {
                    if(j == i)
                    {
                        A[i][j] =
                            -x_[i].boundaryField()[patchI][faceI]
                            / Dij(i,N).boundaryField()[patchI][faceI] / W(N);

                        for(label l=0; l<N+1; l++)
                        {
                            if(l != i)
                            {
                                A[i][j] -= x_[l].boundaryField()[patchI][faceI]
                                    / Dij(i,l).boundaryField()[patchI][faceI]
                                    / W(i);
                            }
                        }

                        B[i][j] = -(x_[i].boundaryField()[patchI][faceI] / W(N)
                            + (1-x_[i].boundaryField()[patchI][faceI]) / W(i));
                    }
                    else
                    {
                        A[i][j] = x_[i].boundaryField()[patchI][faceI]
                            * (1/Dij(i,j).boundaryField()[patchI][faceI] / W(j)
                           - 1/Dij(i,N).boundaryField()[patchI][faceI] / W(N));

                        B[i][j] = x_[i].boundaryField()[patchI][faceI]
                            * (1/W(j) - 1/W(N));
                    }
                }
            }

            scalarRectangularMatrix invA = SVDinv(A);
            multiply(C,invA, B);
            scalarRectangularMatrix invC = SVDinv(C);

            // Prevent spurious numerical oscillations due to the matrix inversion
            for(label i=0; i<N; i++)
            {
                for(label j=0; j<N; j++)
                {
                    D(i,j).boundaryFieldRef()[patchI][faceI] = thermo_.rho()().boundaryField()[patchI][faceI] * C[i][j];
                    G(i,j).boundaryFieldRef()[patchI][faceI] = invC[i][j]/thermo_.rho()().boundaryField()[patchI][faceI];
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::MaxwellStefan<ThermoType>::MaxwellStefan
(
    // A.Alexiou 2014
    // hsCombustionThermo& thermo,
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
    updateMolarFractions();

   D_.setSize((species().size()-1)*(species().size()-1));
   G_.setSize((species().size()-1)*(species().size()-1));

    for(label i=0; i<(species().size()-1); i++)
    {
        for(label j=0; j<(species().size()-1); j++)
        {
            label k = (species().size()-1)*i+j;

            D_.set
            (
                k, new volScalarField
                (
                    IOobject
                    (
                        "D_" + species()[i] + "_" + species()[j],
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar("D", dimensionSet(1, -1, -1, 0, 0), 0.0)
                )
            );

            G_.set
            (
                k, new volScalarField
                (
                    IOobject
                    (
                        "G_" + species()[i] + "_" + species()[j],
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh_,
                    dimensionedScalar("G", dimensionSet(-1, 1, 1, 0, 0), 0.0)
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
Foam::scalar Foam::MaxwellStefan<ThermoType>::correct
(
    fvMatrix<scalar>& yEqn,
    label fieldNumber
)
{
    scalar maxResidual = 0;

    volScalarField& Yi = thermo_.composition().Y(fieldNumber); //yEqn.psi();

    yEqn -= fvm::laplacian(turbulence_.muEff(), Yi);
    yEqn += fvm::laplacian(D(fieldNumber,fieldNumber), Yi, "laplacian(D,Yi)");
    yEqn += Sy_[fieldNumber];

    for(label j=0; j<(species().size()-1); j++)
    {
        if (j != fieldNumber)
        {
            volScalarField& Yj = thermo_.composition().Y(j);

            yEqn += fvc::laplacian(D(fieldNumber,j), Yj, "laplacian(D,Yi)");

        }
    }


    return maxResidual;
}


template<class ThermoType>
bool Foam::MaxwellStefan<ThermoType>::read()
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
