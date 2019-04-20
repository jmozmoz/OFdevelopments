{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) + mvConvection->fvmDiv(phi, he)
      + fvc::ddt(rho, K) + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      - fvm::laplacian(turbulence->alphaEff(), he)
     ==
        Qdot
      - fvc::laplacian(turbulence->alphaEff()*T,
             CpPre)
      + fvOptions(rho, he)
    );

    forAll(Y, i)
    {
        tmp<volScalarField> hsTmp
        (
            new volScalarField(
                IOobject("hsTmp", mesh.time().timeName(), mesh),
                mesh,
                dimensionedScalar("hsTmp", dimEnergy/dimMass, 0)
            )
        );
        volScalarField& hs = hsTmp.ref();

        const volScalarField& p = mesh.lookupObject<volScalarField>("p");

        forAll(hs, cellI)
        {
            hs[cellI] = composition.Hs(i, p[cellI], T[cellI]);
        }

        EEqn -= fvc::laplacian(muEffPre*hs, Y[i]);
    }

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    Info<< "he2: " << he << endl;
    volScalarField tmp(composition.Hs(0, 1e5, 300)*Y[0] + composition.Hs(1, 1e5, 300)*Y[1]);
    Info<< "hx: " << tmp << endl;

    fvOptions.correct(he);

    thermo.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}
