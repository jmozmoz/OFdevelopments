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

Application
    reactingFoam

Description
    Solver for combustion with chemical reactions.

\*---------------------------------------------------------------------------*/
/*Mohsen
#include "fvCFD.H"
#include "hCombustionThermo.H"
#include "turbulenceModel.H"
#include "psiChemistryModel.H"
#include "chemistrySolver.H"
#include "multivariateScheme.H"*/


#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "pressureControl.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"


// ----------------------------- code addition ----------------------------- //
#include "multiSpeciesTransportModel.H"
// ------------------------------------------------------------------------- //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// A.Alexiou 2015, to use OF 2.1 solver style define RHO_EQN_2_1
//#define RHO_EQN_2_1


int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createTimeControls.H"
    #include "initContinuityErrs.H"

///////////////////////////////////////////////////////////////////////////////////////////////////////////
//#   include "readChemistryProperties.H"//////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

#   include "readGravitationalAcceleration.H"

///////////////////////////////////////////////////////////////////////////////////////////////////////////
#   include "createFields.H"/////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

    turbulence->validate();

#   include "initContinuityErrs.H"
#   include "readTimeControls.H"
#   include "compressibleCourantNo.H"
#   include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
#       include "readTimeControls.H"
#       include "compressibleCourantNo.H"
#       include "setDeltaT.H"
        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

///////////////////////////////////////////////////////////////////////////////////////////////////////////
//#       include "chemistry.H"////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////


#include "rhoEqn.H" // OF header/source

        while (pimple.loop()) //Mohsen
        {
            #include "UEqn.H"
            #include "YEqn.H"
            #include "hsEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        // A.Alexiou 2015
        //runTime.write();
        if (runTime.write())
        {
            // A.Alexiou 2015
            thermo.T().write();
//            chemistry.dQ()().write();
            // A.Alexiou 2015
            forAll(composition.Y(), i)
            {
                composition.Y()[i].write();
            }
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
