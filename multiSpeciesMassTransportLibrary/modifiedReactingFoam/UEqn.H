// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(rho, U) + fvm::div(phi, U)
  + MRF.DDt(rho, U)
  + turbulence->divDevRhoReff(U)
 ==
    fvOptions(rho, U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

//    fvVectorMatrix UEqn
//            (
//                fvm::ddt(rho, U)
//                + fvm::div(phi, U)
//                + turbulence->divDevRhoReff(U)
//                ==
//                rho*g
//                );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    /* if (momentumPredictor)  Mohsen
    {
        solve(UEqn == -fvc::grad(p));
    }*/
    if (pimple.momentumPredictor()) //Mohsen
    {
        solve(UEqn == -fvc::grad(p));
     
        fvOptions.correct(U);
        K = 0.5*magSqr(U);
    }
