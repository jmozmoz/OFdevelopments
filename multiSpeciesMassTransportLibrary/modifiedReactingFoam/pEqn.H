
rho = thermo.rho();

volScalarField rAU(1.0/UEqn.A());
//surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));  // OF 2.3

//volVectorField HbyA("HbyA",U);  // OF 2.3
//HbyA = rAU*UEqn.H(); // OF 2.3

surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

if (pimple.transonic())
{

    //               // 2.1
    //               surfaceScalarField phid
    //               (
    //                   "phid",
    //                   fvc::interpolate(psi)
    //                  *(
    //                       (fvc::interpolate(U) & mesh.Sf())
    //                     + fvc::ddtPhiCorr(rAU, rho, U, phi)
    //                   )
    //               );


    // 2.3
//    surfaceScalarField phid
//            (
//                "phid",
//                fvc::interpolate(psi)
//                *(
//                    (fvc::interpolate(rho*HbyA) & mesh.Sf())
//                    + rhorAUf*fvc::ddtCorr(rho, U, phi)
//                    )/fvc::interpolate(rho)
//                );

    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(psi)
       *(
            fvc::flux(HbyA)
          + MRF.zeroFilter
            (
                rhorAUf*fvc::ddtCorr(rho, U, phi)/fvc::interpolate(rho)
            )
        )
    );
    
    MRF.makeRelative(fvc::interpolate(psi), phid);


    //               // 2.1
    //               while (pimple.correctNonOrthogonal())
    //               {
    //                   fvScalarMatrix pEqn
    //                   (
    //                       fvm::ddt(psi, p)
    //                     + fvm::div(phid, p)
    //                     - fvm::laplacian(rho*rAU, p)
    //                   );

    //               pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

    //               if (pimple.finalNonOrthogonalIter())
    //                   {
    //                       phi == pEqn.flux();
    //                   }
    //               }


    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
                (
                    fvm::ddt(psi, p)
                    + fvm::div(phid, p)
                    - fvm::laplacian(rhorAUf, p)
                    ==
                    fvOptions(psi, p, rho.name())
                    );

//        fvOptions.constrain(pEqn);

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi == pEqn.flux();
        }
    }


}
else
{

    //               // 2.1
    //               phi =
    //                   fvc::interpolate(rho)
    //                  *(
    //                       (fvc::interpolate(U) & mesh.Sf())
    //                     + fvc::ddtPhiCorr(rAU, rho, U, phi)
    //                   );


    // 2.3
//    surfaceScalarField phiHbyA
//            (
//                "phiHbyA",
//                (
//                    (fvc::interpolate(rho*HbyA) & mesh.Sf())
//                    + rhorAUf*fvc::ddtCorr(rho, U, phi)
//                    )
//                );
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
            fvc::flux(rho*HbyA)
          + MRF.zeroFilter(rhorAUf*fvc::ddtCorr(rho, U, phi))
        )
    );

    MRF.makeRelative(fvc::interpolate(rho), phiHbyA);


    //               // 2.1
    //               while (pimple.correctNonOrthogonal())
    //               {
    //                   fvScalarMatrix pEqn
    //                   (
    //                       fvm::ddt(psi, p)
    //                     + fvc::div(phi)
    //                     - fvm::laplacian(rho*rAU, p)
    //                   );

    //               pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

    //               if (pimple.finalNonOrthogonalIter())
    //                   {
    //                       phi += pEqn.flux();
    //                   }
    //               }

    // 2.3
    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
                (
                    fvm::ddt(psi, p)
                    + fvc::div(phiHbyA)
          - fvm::laplacian(rhorAUf, p)
                    ==
                    fvOptions(psi, p, rho.name())
                    );

//        fvOptions.constrain(pEqn);

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
    }



}

//#ifndef RHO_EQN_2_1
//
//#include "rhoEqn.H"
//
//#else

// A.Alexiou 2015
// original rhoEqn.H from OF 2.1, works fine with OF 2.3, but then rhoEqn.H has to be
// commented and 'fv::IOoptionList fvOptions(mesh);' is not used
//{
//    solve(fvm::ddt(rho) + fvc::div(phi));
//}

//#endif

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

//           // 2.1
//           #include "compressibleContinuityErrs.H"

//           U -= rAU*fvc::grad(p);
//           U.correctBoundaryConditions();

//           K = 0.5*magSqr(U);

//           //DpDt = fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p);
//           dpdt = fvc::ddt(p);



U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);
K = 0.5*magSqr(U);

if (pressureControl.limit(p))
{
    p.correctBoundaryConditions();
    rho = thermo.rho();
}

if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);
}
