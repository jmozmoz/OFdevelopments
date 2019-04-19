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
      - fvc::laplacian(turbulence->alphaEff(), he)
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
            hs[cellI] = composition.Hs(i, p[cellI], T[cellI]);// * rho[cellI]; // * composition.rho(i, p[cellI], T[cellI]);
        }

        volScalarField::Boundary& hsBf = hs.boundaryFieldRef();

        forAll(hsBf, patchi)
        {
            scalarField& hsp = hsBf[patchi];
            const scalarField& pp = p.boundaryField()[patchi];
            const scalarField& Tup = T.boundaryField()[patchi];

            forAll(hsp, facei)
            {
                hsp[facei] = composition.Hs(i, pp[facei], Tup[facei]);// * rhop[facei]; // composition.rho(i, p[facei], T[facei]);
            }
        }

//        Info << "ha Y" << i << ha << endl;
//        Info << "hs Y" << i << hs << endl;
//        Info << "hx Y" << i << hx << endl;

//        EEqn += fvc::laplacian(turbulence->muEff()*ha, Y[i]);
//        Info<<"Y["<<i<<"].old: "<< Y[i].oldTime() << endl;
//        Info<<"Y["<<i<<"]:     "<< Y[i] << endl;

        volScalarField tmp(hs*(Y[i].oldTime() - Y[i]));
        Info<<"delta h: " << tmp<<endl;
        EEqn += rhoPre*hs*(Y[i].oldTime() - Y[i])/runTime.deltaT();
    }

//    Info<< "muEff: " << turbulence->muEff() << endl;
//    Info<< "muEffPre: " << muEffPre << endl;
    volScalarField he1(he);
//    Info<< "he1: " << he << endl;
//    Info<<"rhoPre: "<<rhoPre<<endl;
//    Info<<"rho:    "<<rho<<endl;
//


    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

//    Info<< "T2: " << T << endl;
    Info<< "he2: " << he << endl;
    Info<< "hex: " << he - he1 << endl;
    volScalarField tmp(composition.Hs(0, 1e5, 300)*Y[0] + composition.Hs(1, 1e5, 300)*Y[1]);
    Info<< "hx: " << tmp << endl;

//    Info<< "Cp * T: " << T * thermo.Cp() << endl;

//    Info<< "diff2 cp: "
//        << fvc::laplacian(turbulence->alphaEff(), he)
//         - fvc::laplacian(turbulence->alphaEff()*T,
//               CpPre) << endl;



    fvOptions.correct(he);

    thermo.correct();

    Info<< "min/max(T) = "
        << min(T).value() << ", " << max(T).value() << endl;
}

