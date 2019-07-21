/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "wallCondensationCoupledMixedFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "rhoReactionThermo.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

volScalarField&
wallCondensationCoupledMixedFvPatchScalarField::outputScalarField
(
    const word& fieldName,
    const dimensionSet& dimSet,
    const fvMesh& mesh
)
{
    if (!mesh.foundObject<volScalarField>(fieldName))
    {
        tmp<volScalarField> tField
        (
            new volScalarField
            (
                IOobject
                (
                    fieldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimSet, Zero)
            )
        );

        tField.ptr()->store();
    }

    return
        const_cast<volScalarField &>
        (
            mesh.lookupObject<volScalarField>(fieldName)
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

wallCondensationCoupledMixedFvPatchScalarField::
wallCondensationCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), "undefined", "undefined", "undefined-K"),
    TnbrName_("undefined-Tnbr"),
    qrNbrName_("undefined-qrNbr"),
    qrName_("undefined-qr"),
    thicknessLayers_(0),
    kappaLayers_(0),
    contactRes_(0),
    specieName_("undefined-specieName"),
    liquid_(nullptr),
    liquidDict_(nullptr),
    mass_(patch().size(), Zero),
    massOld_(patch().size(), Zero),
    myKDelta_(patch().size(), Zero),
    dmHfg_(patch().size(), Zero),
    mpCpTp_(patch().size(), Zero),
    Mcomp_(0.0),
    fluid_(false),
    thickness_(patch().size(), Zero),
    lastTimeStep_(0)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 1.0;
}


wallCondensationCoupledMixedFvPatchScalarField::
wallCondensationCoupledMixedFvPatchScalarField
(
    const wallCondensationCoupledMixedFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(psf, p, iF, mapper),
    temperatureCoupledBase(patch(), psf),
    TnbrName_(psf.TnbrName_),
    qrNbrName_(psf.qrNbrName_),
    qrName_(psf.qrName_),
    thicknessLayers_(psf.thicknessLayers_),
    kappaLayers_(psf.kappaLayers_),
    contactRes_(psf.contactRes_),
    specieName_(psf.specieName_),
    liquid_(psf.liquid_),
    liquidDict_(psf.liquidDict_),
    mass_(psf.mass_, mapper),
    massOld_(psf.massOld_, mapper),
    myKDelta_(psf.myKDelta_, mapper),
    dmHfg_(psf.dmHfg_, mapper),
    mpCpTp_(psf.mpCpTp_, mapper),
    Mcomp_(psf.Mcomp_),
    fluid_(psf.fluid_),
    thickness_(psf.thickness_, mapper),
    lastTimeStep_(psf.lastTimeStep_)
{}


wallCondensationCoupledMixedFvPatchScalarField::
wallCondensationCoupledMixedFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    temperatureCoupledBase(patch(), dict),
    TnbrName_(dict.lookupOrDefault<word>("Tnbr", "T")),
    qrNbrName_(dict.lookupOrDefault<word>("qrNbr", "none")),
    qrName_(dict.lookupOrDefault<word>("qr", "none")),
    thicknessLayers_(0),
    kappaLayers_(0),
    contactRes_(0.0),
    specieName_(dict.lookupOrDefault<word>("specie", "undefined-specieName")),
    liquid_(nullptr),
    liquidDict_(),
    mass_(patch().size(), Zero),
    massOld_(patch().size(), Zero),
    myKDelta_(patch().size(), Zero),
    dmHfg_(patch().size(), Zero),
    mpCpTp_(patch().size(), Zero),
    Mcomp_(dict.lookupOrDefault<scalar>("carrierMolWeight", 0.0)),
    fluid_(false),
    thickness_(patch().size(), Zero),
    lastTimeStep_(0){
    if (!isA<mappedPatchBase>(this->patch().patch()))
    {
        FatalErrorInFunction
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.name()
            << " of field " << internalField().name()
            << " in file " << internalField().objectPath()
            << exit(FatalError);
    }

    if (dict.found("thicknessLayers"))
     {
        dict.lookup("thicknessLayers") >> thicknessLayers_;
        dict.lookup("kappaLayers") >> kappaLayers_;

        if (thicknessLayers_.size() > 0)
        {
            // Calculate effective thermal resistance by harmonic averaging
            forAll(thicknessLayers_, iLayer)
            {
                contactRes_ += thicknessLayers_[iLayer]/kappaLayers_[iLayer];
            }
            contactRes_ = 1.0/contactRes_;
        }
    }

    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("refValue"))
    {
        // Full restart
        refValue() = scalarField("refValue", dict, p.size());
        refGrad() = scalarField("refGradient", dict, p.size());
        valueFraction() = scalarField("valueFraction", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }

    if (dict.found("specie"))
    {
        fluid_ = true;
        Info << "************************ found fluid side" << endl;
    }

    if (fluid_)
    {
        liquidDict_ = dict.subDict("liquid");
        liquid_ =
            liquidProperties::New(liquidDict_.subDict(specieName_));

        if (dict.found("thickness"))
        {
            scalarField& Tp = *this;
            const scalarField& magSf = patch().magSf();

            // Assume initially standard pressure for rho calculation
            scalar pf = 1e5;
            thickness_ = scalarField("thickness", dict, p.size());
            forAll(thickness_, i)
            {
                mass_[i] =
                    thickness_[i]*liquid_->rho(pf, Tp[i])*magSf[i];
                massOld_[i] = mass_[i];
            }
        }
        lastTimeStep_ = patch().boundaryMesh().mesh().time().value();
    }
}


wallCondensationCoupledMixedFvPatchScalarField::
wallCondensationCoupledMixedFvPatchScalarField
(
    const wallCondensationCoupledMixedFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(psf, iF),
    temperatureCoupledBase(patch(), psf),
    TnbrName_(psf.TnbrName_),
    qrNbrName_(psf.qrNbrName_),
    qrName_(psf.qrName_),
    thicknessLayers_(psf.thicknessLayers_),
    kappaLayers_(psf.kappaLayers_),
    contactRes_(psf.contactRes_),
    specieName_(psf.specieName_),
    liquid_(psf.liquid_, false), // do not re-use it but use clone internally
    liquidDict_(psf.liquidDict_),
    mass_(psf.mass_),
    massOld_(psf.massOld_),
    myKDelta_(psf.myKDelta_),
    dmHfg_(psf.dmHfg_),
    mpCpTp_(psf.mpCpTp_),
    Mcomp_(psf.Mcomp_),
    fluid_(psf.fluid_),
    thickness_(psf.thickness_),
    lastTimeStep_(psf.lastTimeStep_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void wallCondensationCoupledMixedFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);

    if (fluid_)
    {
        mass_.autoMap(m);
        massOld_.autoMap(m);
        myKDelta_.autoMap(m);
        dmHfg_.autoMap(m);
        mpCpTp_.autoMap(m);
        thickness_.autoMap(m);
    }
}


void wallCondensationCoupledMixedFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const wallCondensationCoupledMixedFvPatchScalarField& tiptf =
        refCast<const wallCondensationCoupledMixedFvPatchScalarField>
        (
            ptf
        );

    if (fluid_)
    {
        mass_.rmap(tiptf.mass_, addr);
        massOld_.rmap(tiptf.mass_, addr);
        myKDelta_.rmap(tiptf.myKDelta_, addr);
        dmHfg_.rmap(tiptf.dmHfg_, addr);
        mpCpTp_.rmap(tiptf.mpCpTp_, addr);
        thickness_.rmap(tiptf.thickness_, addr);
    }
}

void wallCondensationCoupledMixedFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(patch().patch());

    const scalarField& magSf = patch().magSf();

    const polyMesh& mesh = patch().boundaryMesh().mesh();

    const polyMesh& nbrMesh = mpp.sampleMesh();
    const label samplePatchi = mpp.samplePolyPatch().index();
    const fvPatch& nbrPatch =
        refCast<const fvMesh>(nbrMesh).boundary()[samplePatchi];


    scalarField Tinternal(patchInternalField());
    scalarField& Tpatch = *this;

    const volScalarField& TmeshField = static_cast<const volScalarField&>
    (
        internalField()
    );
    scalarField TpatchOld
    (
        TmeshField.oldTime().boundaryField()[patch().index()]
    );

    typedef wallCondensationCoupledMixedFvPatchScalarField thisType;

    const fvPatchScalarField& nbrTp =
        nbrPatch.lookupPatchField<volScalarField, scalar>(TnbrName_);

    if (!isA<thisType>(nbrTp))
    {
        FatalErrorInFunction
            << "Patch field for " << internalField().name() << " on "
            << patch().name() << " is of type " << thisType::typeName
            << endl << "The neighbouring patch field " << TnbrName_ << " on "
            << nbrPatch.name() << " is required to be the same, but is "
            << "currently of type " << nbrTp.type() << exit(FatalError);
    }

    const thisType& nbrField = refCast<const thisType>(nbrTp);

    // Swap to obtain full local values of neighbour internal field
    scalarField nbrInternalField(nbrField.patchInternalField());
    mpp.distribute(nbrInternalField);

    myKDelta_ = kappa(Tpatch)*patch().deltaCoeffs();

    scalarField dm(patch().size(), Zero);
    scalarField dHspec(patch().size(), Zero);
    scalarField hPhaseChange(patch().size(), Zero);
    scalarField dmhPhaseChange(patch().size(), Zero);
    scalarField hRemovedMass(patch().size(), Zero);
    scalarField dmhRemoveMass(patch().size(), Zero);

    // Fluid Side
    if (fluid_)
    {
        const scalar dt = mesh.time().deltaTValue(); // @suppress("Invalid arguments") because I know better

        if(mesh.time().value() != lastTimeStep_)
        {
            lastTimeStep_= mesh.time().value();
            massOld_ = mass_;
        }

        const rhoReactionThermo& thermo =
            mesh.lookupObject<rhoReactionThermo>
            (
                basicThermo::dictName
            );

        const basicSpecieMixture& composition = thermo.composition();
        label specieIndex = composition.species()[specieName_];

        volScalarField & filmMassSource =
            const_cast<volScalarField &>
            (
                mesh.lookupObject<volScalarField>
                (
                    "filmMassSource" + specieName_
                )
            );

        volScalarField & filmEnergySource =
            const_cast<volScalarField &>
            (
                mesh.lookupObject<volScalarField>
                (
                    "filmEnergySource" + specieName_
                )
            );

        const scalarField myDelta(patch().deltaCoeffs());

        scalarField cp(patch().size(), Zero);
        scalarField liquidRho(patch().size(), Zero);

        const fvPatchScalarField& Ypatch =
            patch().lookupPatchField<volScalarField, scalar>
            (
                specieName_
            );

        const fvPatchScalarField& pPatch =
            patch().lookupPatchField<volScalarField, scalar>("p");

        const fvPatchScalarField& rhoPatch =
            patch().lookupPatchField<volScalarField, scalar>("rho");

        const fvPatchScalarField& muPatch =
            patch().lookupPatchField<volScalarField, scalar>("thermo:mu");

        const fvPatchScalarField& nutPatch =
            patch().lookupPatchField<volScalarField, scalar>("nut");

        const scalarField Yinternal(Ypatch.patchInternalField());
        const scalarField pInternal(pPatch.patchInternalField());
        const scalarField rhoInternal(rhoPatch.patchInternalField());
        const labelList& faceCells = patch().faceCells();

        const scalar Sct = 0.9;
        const scalar Sc = 1.0;

        forAll(Tpatch, faceI)
        {
            const scalar Tface = Tpatch[faceI];
            const scalar Tcell = Tinternal[faceI];
            const scalar pFace = pPatch[faceI];
            const scalar pCell = pInternal[faceI];

            const scalar muFace = muPatch[faceI];
            const scalar nuCell = 3.45e-5;
            const scalar rhoFace = rhoPatch[faceI];
            const scalar rhoCell = rhoInternal[faceI];
            const scalar muCell = nuCell * rhoCell;
            const scalar mutFace = nutPatch[faceI]*rhoFace;
            const scalar pSatCell = liquid_->pv(pFace, Tcell);
            const scalar pSatFace = liquid_->pv(pFace, Tface);
            const scalar Mv = liquid_->W();
            const scalar Ycell = Yinternal[faceI];
            const scalar deltaFace = myDelta[faceI];

            cp[faceI] = liquid_->Cp(pFace, Tface);
            hPhaseChange[faceI] = liquid_->hl(pFace, Tface);
            hRemovedMass[faceI] = composition.Hs(specieIndex, pCell, Tcell);

            // Calculate relative humidity
            const scalar invMwmean =
                Yinternal[faceI]/Mv + (1.0 - Ycell)/Mcomp_;
            const scalar Xv = Ycell/invMwmean/Mv;
            const scalar RH = min(Xv*pFace/pSatCell, 1.0);

            scalar RHmin = 0.01;
            scalar Tdew = -GREAT;

            if (RH > RHmin)
            {
                scalar b = 243.5;
                scalar c = 17.65;
                scalar TintDeg = Tcell - 273;
                Tdew =
                    b*(log(RH) + (c*TintDeg)/(b + TintDeg))
                   /(c - log(RH) - ((c*TintDeg)/(b + TintDeg))) + 273;
            }

            if
            (
                Tface < Tdew
             && RH > RHmin
            )
            {
                const scalar YsatFace = pSatFace/pFace*Mv/Mcomp_;

                const scalar gamma = muCell + mutFace / Sct;

                // mass flux [kg/s/m^2]
                // positive if mass is condensing
                dm[faceI] =
                    gamma*deltaFace
                   *(Ycell - YsatFace)/(1 - YsatFace);
            }

            liquidRho[faceI] = liquid_->rho(pFace, Tface);
        }

        // Output film delta (e.g. H2OThickness) [m]
        scalarField &thicknessField =
            outputScalarField
            (
                specieName_ + "Thickness",
                dimLength,
                refCast<const fvMesh>(mesh)
            ).boundaryFieldRef()[patch().index()];
            thicknessField = mass_/liquidRho/magSf;

        scalarField &massFluxOut =
            outputScalarField
            (
                specieName_ + "MassFlux",
                dimMass/dimArea/dimTime,
                refCast<const fvMesh>(mesh)
            ).boundaryFieldRef()[patch().index()];
        massFluxOut = dm;

//            scalarField &heatFluxOut =
//                 outputScalarField
//                 (
//                     specieName_ + "MassFlux",
//                     dimMass/dimArea/dimTime,
//                     refCast<const fvMesh>(mesh)
//                 ).boundaryFieldRef()[patch().index()];
//            heatFluxOut = dm;

        mass_= massOld_ + dm*dt*magSf;
        mass_ = max(mass_, scalar(0));
        mpCpTp_ = mass_*cp/dt/magSf;

        // Heat flux due to change of phase [W/m2]
        dmHfg_ = dm*hPhaseChange;
        dHspec = dm*hRemovedMass;

        forAll(faceCells, faceI)
        {
            const label cellI = faceCells[faceI];

            filmMassSource[cellI] =
               -dm[faceI]
               *magSf[faceI]
               /mesh.cellVolumes()[cellI];

            filmEnergySource[cellI] =
               -dm[faceI]*hRemovedMass[faceI]
               *magSf[faceI]
               /mesh.cellVolumes()[cellI];
        }
    }

    // Swap to obtain full local values of neighbour K*delta
    scalarField KDeltaNbr;
    if (contactRes_ == 0.0)
    {
        KDeltaNbr = nbrField.kappa(nbrField)*nbrPatch.deltaCoeffs();
    }
    else
    {
        KDeltaNbr.setSize(nbrField.size(), contactRes_);
    }
    mpp.distribute(KDeltaNbr);

    scalarField mpCpTpNbr(patch().size(), Zero);
    scalarField dmHfgNbr(patch().size(), Zero);

    if (!fluid_)
    {
        mpCpTpNbr = nbrField.mpCpTp();
        mpp.distribute(mpCpTpNbr);

        dmHfgNbr = nbrField.dmHfg();
        mpp.distribute(dmHfgNbr);
    }

    scalarField qr(Tpatch.size(), 0.0);
    if (qrName_ != "none")
    {
        qr = patch().lookupPatchField<volScalarField, scalar>(qrName_);
    }

    scalarField qrNbr(Tpatch.size(), 0.0);
    if (qrNbrName_ != "none")
    {
        qrNbr = nbrPatch.lookupPatchField<volScalarField, scalar>(qrNbrName_);
        mpp.distribute(qrNbr);
    }

    const scalarField dmHfg(dmHfgNbr + dmHfg_);
    const scalarField mpCpdt(mpCpTpNbr + mpCpTp_);

    // qr > 0 (heat up the wall)
    scalarField alpha(KDeltaNbr + mpCpdt - (qr + qrNbr)/Tpatch);

    valueFraction() = alpha/(alpha + myKDelta_);
    refValue() = (KDeltaNbr*nbrInternalField + mpCpdt*TpatchOld + dmHfg)/alpha;

    mixedFvPatchScalarField::updateCoeffs();

    if (fluid_)
    {
        scalar Qdm = gSum(dm*magSf);
        scalar QMass = gSum(mass_);
        scalar Qt = gSum(myKDelta_*(Tpatch - Tinternal)*magSf);
        scalar QtSolid = gSum(KDeltaNbr*(Tpatch - nbrInternalField)*magSf);

        scalar Q = gSum(kappa(Tpatch)*patch().magSf()*snGrad());
        scalar Qhgf = gSum(dmHfg_*magSf);
        scalar Qspec = gSum(dHspec*magSf);

        Info<< mesh.name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " <- "
            << nbrMesh.name() << ':'
            << nbrPatch.name() << ':'
            << this->internalField().name() << " :" << nl
            << " heat transfer rate:" << Q << endl
            << "    Total mass flux [Kg/s]:                            "
            << Qdm << nl
            << "    Total mass on the wall [Kg]:                       "
            << QMass << nl
            << "    Total heat (>0 leaving the wall to the fluid) [W]: "
            << Qt << nl
            << "    Total heat (>0 leaving the wall to the solid) [W]: "
            << QtSolid << nl
            << "    Total latent heat released [W]:                    "
            << Qhgf << nl
            << "    Total specific heat removed from fluid [W]:        "
            << Qspec << nl
            << " walltemperature "
            << " min:" << gMin(Tpatch)
            << " max:" << gMax(Tpatch)
            << " avg:" << gAverage(Tpatch)
            << endl;
    }

    // Restore tag
    UPstream::msgType() = oldTag;
}


void wallCondensationCoupledMixedFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    os.writeKeyword("Tnbr")<< TnbrName_ << token::END_STATEMENT << nl;
    os.writeKeyword("qrNbr")<< qrNbrName_ << token::END_STATEMENT << nl;
    os.writeKeyword("qr")<< qrName_ << token::END_STATEMENT << nl;
    thicknessLayers_.writeEntry("thicknessLayers", os);
    kappaLayers_.writeEntry("kappaLayers", os);

    if (fluid_)
    {
        os.writeKeyword("specie") << specieName_ << token::END_STATEMENT << nl;
        os.writeKeyword("carrierMolWeight") << Mcomp_ << token::END_STATEMENT << nl;
        mass_.writeEntry("mass", os);
        thickness_.writeEntry("thickness", os);
        word liq = "liquid";
        os << token::TAB << token::TAB << liq;
        liquidDict_.write(os);
    }

    temperatureCoupledBase::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    wallCondensationCoupledMixedFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam


// ************************************************************************* //
