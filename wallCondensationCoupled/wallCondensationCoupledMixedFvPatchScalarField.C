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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

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
    Mcomp_(0.0),
    fluid_(false),
    thickness_(patch().size(), Zero)
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
    Mcomp_(psf.Mcomp_),
    fluid_(psf.fluid_),
    thickness_(psf.thickness_, mapper)
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
    Mcomp_(dict.lookupOrDefault<scalar>("carrierMolWeight", 0.0)),
    fluid_(false),
    thickness_(patch().size(), Zero)
{
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
    liquid_(psf.liquid_),
    liquidDict_(psf.liquidDict_),
    mass_(psf.mass_),
    massOld_(psf.massOld_),
    Mcomp_(psf.Mcomp_),
    fluid_(psf.fluid_),
    thickness_(psf.thickness_)
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

    // Swap to obtain full local values of neighbour internal field
    scalarField nbrIntFld(nbrField.patchInternalField());
    mpp.distribute(nbrIntFld);

    scalarField Tc(patchInternalField());
    scalarField& Tp = *this;

    scalarField Tin(patchInternalField());

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
    scalarField TcNbr(nbrField.patchInternalField());
    mpp.distribute(TcNbr);

    scalarField dm(patch().size(), Zero);

    // Fluid Side
    if (fluid_)
    {
        scalarField Yvp(patch().size(), Zero);
        const scalar dt = mesh.time().deltaTValue();

        const scalarField myDelta(patch().deltaCoeffs());

        scalarField cp(patch().size(), Zero);
        scalarField hfg(patch().size(), Zero);
        scalarField htc(patch().size(), GREAT);
        scalarField liquidRho(patch().size(), Zero);

        const fvPatchScalarField& Yp =
            patch().lookupPatchField<volScalarField, scalar>
            (
                specieName_
            );

        const fvPatchScalarField& pp =
            patch().lookupPatchField<volScalarField, scalar>("p");

        const fvPatchScalarField& rhop =
            patch().lookupPatchField<volScalarField, scalar>("rho");

        const fvPatchScalarField& mup =
            patch().lookupPatchField<volScalarField, scalar>("thermo:mu");

        const scalarField Yi(Yp.patchInternalField());

        forAll(Tp, faceI)
        {
            const scalar Tf = Tp[faceI];
            const scalar Tint = Tin[faceI];
            const scalar pf = pp[faceI];

            const scalar muf = mup[faceI];
            const scalar rhof = rhop[faceI];
            const scalar nuf = muf/rhof;
            const scalar pSat = liquid_->pv(pf, Tint);
            const scalar Mv = liquid_->W();
            const scalar TSat = liquid_->pvInvert(pSat);

            cp[faceI] = liquid_->Cp(pf, Tf);
            hfg[faceI] = liquid_->hl(pf, Tf);

            // Calculate relative humidity
            const scalar invMwmean =
                    Yi[faceI]/Mv + (1.0 - Yi[faceI])/Mcomp_;
            const scalar Xv = Yi[faceI]/invMwmean/Mv;
            const scalar RH = min(Xv*pf/pSat, 1.0);

            scalar RHmin = 0.01;
            scalar Tdew = -GREAT;

            if (RH > RHmin)
            {
                scalar b = 243.5;
                scalar c = 17.65;
                scalar TintDeg = Tint - 273;
                Tdew =
                    b*(log(RH) + (c*TintDeg)/(b + TintDeg))
                   /(c - log(RH) - ((c*TintDeg)/(b + TintDeg))) + 273;
            }

            if
            (
                Tf < Tdew
             && RH > RHmin
            )
            {
//                    htc[faceI] = htcCondensation(TSat, Re)*nbrK[faceI]/L_;

//                    scalar htcTotal =
//                        1.0/((1.0/myKDelta_[faceI]) + (1.0/htc[faceI]));

                // Heat flux [W] (>0 heat is converted into mass)
//                    const scalar q = (Tint - Tf)*htcTotal*magSf[faceI];

                // Mass flux rate [Kg/s/m2]
//                    dm[faceI] = q/hfg[faceI]/magSf[faceI];

//                    mass_[faceI] += q/hfg[faceI]*dt;

                // -dYp/dn = q/Dab (fixedGradient)
                const scalar Dab = liquid_->D(pf, Tf);
//                    Yvp[faceI] =
//                        -min(dm[faceI]/Dab/rhof, Yi[faceI]*myDelta[faceI]);
            }

            liquidRho[faceI] = liquid_->rho(pf, Tf);

            mass_ = max(mass_, scalar(0));

            // Output film delta (e.g. H2OThickness) [m]
            const word fieldName(specieName_ + "Thickness");

//            scalarField& pDelta =
//                thicknessField
//                (
//                    fieldName,
//                    refCast<const fvMesh>(mesh)
//                ).boundaryFieldRef()[patch().index()];
//
//
//            pDelta = mass_/liquidRho/magSf;
//
//            // Weight myKDelta and htc
//            myKDelta_ = 1.0/((1.0/myKDelta_) + (1.0/htc));
//
//            mpCpTp_ = mass_*cp/dt/magSf;
//
//            // Heat flux due to change of phase [W/m2]
//            dmHfg_ = dm*hfg;
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

    scalarField KDelta(kappa(Tp)*patch().deltaCoeffs());

    scalarField qr(Tp.size(), 0.0);
    if (qrName_ != "none")
    {
        qr = patch().lookupPatchField<volScalarField, scalar>(qrName_);
    }

    scalarField qrNbr(Tp.size(), 0.0);
    if (qrNbrName_ != "none")
    {
        qrNbr = nbrPatch.lookupPatchField<volScalarField, scalar>(qrNbrName_);
        mpp.distribute(qrNbr);
    }

    valueFraction() = KDeltaNbr/(KDeltaNbr + KDelta);
    refValue() = TcNbr;
    refGrad() = (qr + qrNbr)/kappa(Tp);

    mixedFvPatchScalarField::updateCoeffs();

    if fluid_)
    {
        scalar Qdm = gSum(dm);
        scalar QMass = gSum(mass_);
        scalar Qt = gSum(myKDelta_*(Tp - Tin)*magSf);
        scalar QtSolid = gSum(KDeltaNbr*(Tp - nbrIntFld)*magSf);
        
        scalar Q = gSum(kappa(Tp)*patch().magSf()*snGrad());

        Info<< mesh.name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " <- "
            << nbrMesh.name() << ':'
            << nbrPatch.name() << ':'
            << this->internalField().name() << " :"
            << " heat transfer rate:" << Q
            << "    Total mass flux   [Kg/s] : " << Qdm << nl
            << "    Total mass on the wall [Kg] : " << QMass << nl
            << "    Total heat (>0 leaving the wall to the fluid) [W] : "
            << Qt << nl
            << "     Total heat (>0 leaving the wall to the solid) [W] : "
            << QtSolid << nl
            << " walltemperature "
            << " min:" << gMin(Tp)
            << " max:" << gMax(Tp)
            << " avg:" << gAverage(Tp)
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
