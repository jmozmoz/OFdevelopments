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

    /* Corrections according to presentation of V. Kumar G et al.
     * Development and verification of a multi-species gas transport solver
     * OpenFOAM Workshop 2019, Duisburg, Germany
     */

    // correct -fvm::laplacian(turbulence->alphaEff(), he)
    // = -\nabla (\lapha_{eff} \nabla h)
    // to only include diffusive heat conduction due to temperature
    // gradient: -\nabla (\lambda_{eff} \nabla T)
    if (thermoTypeThermo == "hConst")
    {
        // if Cp is temperature independent, the product rule can be applied
        // to \nabla he and the correction is>
        EEqn += fvc::laplacian(turbulence->alphaEff()*T, thermo.Cp());
    }
    else
    {
        // if Cp is temperature dpeendent, it is necessary to replace
        // laplacian alphaEff * he with laplacian lambdaEff T in the
        // energy equation
        EEqn -= fvc::laplacian(turbulence->alphaEff()*thermo.Cp(), thermo.T());
        EEqn += fvc::laplacian(turbulence->alphaEff(), he);
    }

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

        // the face values must be set, too. Otherwise diffusion
        // will be calcualted correctly, if fixed values are set
        // on boundaries
        volScalarField::Boundary& hsBf = hs.boundaryFieldRef();

        forAll(hsBf, patchi)
        {
            scalarField& hsPatch = hsBf[patchi];
            const scalarField& pp = p.boundaryField()[patchi];
            const scalarField& Tp = T.boundaryField()[patchi];

            forAll(hsPatch, facei)
            {
                hsPatch[facei] = composition.Hs(i, pp[facei], Tp[facei]);
            }
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
