{
#if 0
    fvScalarMatrix hsEqn
            (
                fvm::ddt(rho, hs)
                + mvConvection->fvmDiv(phi, hs)
                - fvm::laplacian(turbulence->alphaEff(), hs)
                + mstm().multiSpeciesHeatSource()
                ==
                dpdt
                // ------------------------------------------------------------------------- //
                + chemistrySh
                // ------------------------------------------------------------------------- //
                );

    hsEqn.relax();

    fvOptions.constrain(hsEqn);

    hsEqn.solve();
#endif
    fvOptions.correct(hs);

//    thermo.correct();
    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}
