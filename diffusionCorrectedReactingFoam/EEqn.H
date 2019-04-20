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
      + fvOptions(rho, he)
    );

    /* Corrections according to presentation of S. Kelm et al.
     * Open Source CFD Modeling of Containment Transport and Mixing Processes
     * Development and First Validation of ‘containmentFOAM’
     * International CFD Network Meeting, GRS Garching, March 19th 2019
     */

    // correct -fvm::laplacian(turbulence->alphaEff(), he)
    // = -\nabla (\lapha_{eff} \nabla h)
    // to only include diffusive heat conduction due to temperature
    // gradient: -\nabla (\lambda_{eff} \nabla T)
    EEqn += fvc::laplacian(turbulence->alphaEff()*T, thermo.Cp());

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

        // correct energy to be consistent with species diffusion
        EEqn -= fvc::laplacian(turbulence->muEff()*hs, Y[i]);
    }

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}
