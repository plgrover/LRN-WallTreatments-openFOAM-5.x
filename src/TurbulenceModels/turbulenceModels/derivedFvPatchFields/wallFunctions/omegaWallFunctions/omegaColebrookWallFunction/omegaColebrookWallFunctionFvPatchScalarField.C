/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "omegaColebrookWallFunctionFvPatchScalarField.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "fvMatrix.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar omegaColebrookWallFunctionFvPatchScalarField::tolerance_ = 1e-5;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void omegaColebrookWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch())) {
        FatalErrorInFunction
                << "Invalid wall function specification" << nl
                << "    Patch type for patch " << patch().name()
                << " must be wall" << nl
                << "    Current patch type is " << patch().type() << nl << endl
                << abort(FatalError);
    }
}


void omegaColebrookWallFunctionFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
    os.writeKeyword("beta1") << beta1_ << token::END_STATEMENT << nl;
    os.writeKeyword("blended") << blended_ << token::END_STATEMENT << nl;
}


void omegaColebrookWallFunctionFvPatchScalarField::setMaster()
{
    if (master_ != -1) {
        return;
    }

    const volScalarField& omega =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    label master = -1;
    forAll(bf, patchi) {
        if (isA<omegaColebrookWallFunctionFvPatchScalarField>(bf[patchi])) {
            omegaColebrookWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            if (master == -1) {
                master = patchi;
            }

            opf.master() = master;
        }
    }
}


void omegaColebrookWallFunctionFvPatchScalarField::createAveragingWeights()
{
    const volScalarField& omega =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    const fvMesh& mesh = omega.mesh();

    if (initialised_ && !mesh.changing()) {
        return;
    }

    volScalarField weights
    (
        IOobject
        (
            "weights",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false // do not register
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    DynamicList<label> omegaPatches(bf.size());
    forAll(bf, patchi) {
        if (isA<omegaColebrookWallFunctionFvPatchScalarField>(bf[patchi])) {
            omegaPatches.append(patchi);

            const labelUList& faceCells = bf[patchi].patch().faceCells();
            forAll(faceCells, i) {
                label celli = faceCells[i];
                weights[celli]++;
            }
        }
    }

    cornerWeights_.setSize(bf.size());
    forAll(omegaPatches, i) {
        label patchi = omegaPatches[i];
        const fvPatchScalarField& wf = weights.boundaryField()[patchi];
        cornerWeights_[patchi] = 1.0/wf.patchInternalField();
    }

    G_.setSize(internalField().size(), 0.0);
    omega_.setSize(internalField().size(), 0.0);

    initialised_ = true;
}


omegaColebrookWallFunctionFvPatchScalarField&
omegaColebrookWallFunctionFvPatchScalarField::omegaPatch(const label patchi)
{
    const volScalarField& omega =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    const omegaColebrookWallFunctionFvPatchScalarField& opf =
        refCast<const omegaColebrookWallFunctionFvPatchScalarField>(bf[patchi]);

    return const_cast<omegaColebrookWallFunctionFvPatchScalarField&>(opf);
}


void omegaColebrookWallFunctionFvPatchScalarField::calculateTurbulenceFields
(
    const turbulenceModel& turbModel,
    scalarField& G0,
    scalarField& omega0
)
{
    // accumulate all of the G and omega contributions
    Info<< "calculateTurbulenceFields" << nl << endl;
    forAll(cornerWeights_, patchi) {
        if (!cornerWeights_[patchi].empty()) {
            omegaColebrookWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            const List<scalar>& w = cornerWeights_[patchi];

            opf.calculate(turbModel, w, opf.patch(), G0, omega0);
        }
    }

    // apply zero-gradient condition for omega
    forAll(cornerWeights_, patchi) {
        if (!cornerWeights_[patchi].empty()) {
            omegaColebrookWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);
            opf == scalarField(omega0, opf.patch().faceCells());
        }
    }
}


void omegaColebrookWallFunctionFvPatchScalarField::calculate
(
    const turbulenceModel& turbModel,
    const List<scalar>& cornerWeights,
    const fvPatch& patch,
    scalarField& G0,
    scalarField& omega0
)
{
    const label patchi = patch.index();

    const scalarField& y = turbModel.y()[patchi];

    const scalar Cmu25 = pow025(Cmu_);

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    // Get the kinematic viscosity
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    // Get the eddy visocity
    const tmp<scalarField> tnutw = turbModel.nut(patchi);
    const scalarField& nutw = tnutw();

    // Get the gradient at the wall
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magGradUw(mag(Uw.snGrad()));

    // Set omega and G
    forAll(nutw, facei) {
        const label celli = patch.faceCells()[facei];

        const scalar yPlus = Cmu25*y[facei]*sqrt(k[celli])/nuw[facei];
        const scalar w = cornerWeights[facei];

        // Calculate the effective viscosity and shear velocity
        scalar nueff = nutw[facei] + nuw[facei];
        scalar ustar = sqrt(nueff*magGradUw[facei]);

        // Calculate kr+ - the non-dimensional roughness height
        // EQ (1) in Aupoix 2014
        scalar krPlus = kr_*ustar/nuw[facei];

        if (krPlus < 1.5) {
            // Smooth wall case
            const scalar omegaVis = 6*nuw[facei]/(beta1_*sqr(y[facei]));
            const scalar omegaLog = sqrt(k[celli])/(Cmu25*kappa_*y[facei]);

            // Switching between the laminar sub-layer and the log-region rather
            // than blending has been found to provide more accurate results 
            // over a range of near-wall y+.
            //
            // For backward-compatibility the blending method is provided as an
            // option

            if (blended_) {
                omega0[celli] += w*sqrt(sqr(omegaVis) + sqr(omegaLog));
            }

            if (yPlus > yPlusLam_) {
                if (!blended_) {
                    omega0[celli] += w*omegaLog;
                }

                G0[celli] +=
                    w
                    *(nutw[facei] + nuw[facei])
                    *magGradUw[facei]
                    *Cmu25*sqrt(k[celli])
                    /(kappa_*y[facei]);
            } else {
                if (!blended_) {
                    omega0[celli] += w*omegaVis;
                }
            }
        } else {
            // Rough wall case
            // Calculate the omega in wall units Eq. 25 in Aupoix (2014)
            scalar omegaPlus = (300.0/sqr(krPlus))/tanh(15.0/(4.0*krPlus))
                + (191./krPlus)*(1 - exp(-krPlus/250.));
            omega0[celli] = omegaPlus * sqr(ustar) / nuw[facei];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

omegaColebrookWallFunctionFvPatchScalarField::omegaColebrookWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
    :
    fixedValueFvPatchField<scalar>(p, iF),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8),
    beta1_(0.075),
    blended_(false),
    yPlusLam_(nutWallFunctionFvPatchScalarField::yPlusLam(kappa_, E_)),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_(),
    Betak_(0.09),
    kr_(0.002)
{
    checkType();
}


omegaColebrookWallFunctionFvPatchScalarField::omegaColebrookWallFunctionFvPatchScalarField
(
    const omegaColebrookWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
    :
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_),
    beta1_(ptf.beta1_),
    blended_(ptf.blended_),
    yPlusLam_(ptf.yPlusLam_),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_(),
    Betak_(ptf.Betak_),
    kr_(ptf.kr_)
{
    checkType();
}


omegaColebrookWallFunctionFvPatchScalarField::omegaColebrookWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
    :
    fixedValueFvPatchField<scalar>(p, iF, dict),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8)),
    beta1_(dict.lookupOrDefault<scalar>("beta1", 0.075)),
    blended_(dict.lookupOrDefault<Switch>("blended", false)),
    yPlusLam_(nutWallFunctionFvPatchScalarField::yPlusLam(kappa_, E_)),
    Betak_(dict.lookupOrDefault<scalar>("Betak", 0.09)),
    kr_(dict.lookupOrDefault<scalar>("kr",0.002)),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();

    // apply zero-gradient condition on start-up
    this->operator==(patchInternalField());
}


omegaColebrookWallFunctionFvPatchScalarField::omegaColebrookWallFunctionFvPatchScalarField
(
    const omegaColebrookWallFunctionFvPatchScalarField& owfpsf
)
    :
    fixedValueFvPatchField<scalar>(owfpsf),
    Cmu_(owfpsf.Cmu_),
    kappa_(owfpsf.kappa_),
    E_(owfpsf.E_),
    beta1_(owfpsf.beta1_),
    blended_(owfpsf.blended_),
    yPlusLam_(owfpsf.yPlusLam_),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_(),
    Betak_(owfpsf.Betak_),
    kr_(owfpsf.kr_)
{
    checkType();
}


omegaColebrookWallFunctionFvPatchScalarField::omegaColebrookWallFunctionFvPatchScalarField
(
    const omegaColebrookWallFunctionFvPatchScalarField& owfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
    :
    fixedValueFvPatchField<scalar>(owfpsf, iF),
    Cmu_(owfpsf.Cmu_),
    kappa_(owfpsf.kappa_),
    E_(owfpsf.E_),
    beta1_(owfpsf.beta1_),
    blended_(owfpsf.blended_),
    yPlusLam_(owfpsf.yPlusLam_),
    Betak_(owfpsf.Betak_),
    kr_(owfpsf.kr_),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalarField& omegaColebrookWallFunctionFvPatchScalarField::G(bool init)
{
    if (patch().index() == master_) {
        if (init) {
            G_ = 0.0;
        }

        return G_;
    }

    return omegaPatch(master_).G();
}


scalarField& omegaColebrookWallFunctionFvPatchScalarField::omega(bool init)
{
    if (patch().index() == master_) {
        if (init) {
            omega_ = 0.0;
        }

        return omega_;
    }

    return omegaPatch(master_).omega(init);
}


void omegaColebrookWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated()) {
        return;
    }
    Info<< "updateCoeffs" << nl << endl;

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
                                       (
                                           IOobject::groupName
                                           (
                                                   turbulenceModel::propertiesName,
                                                   internalField().group()
                                           )
                                       );

    setMaster();

    if (patch().index() == master_) {
        createAveragingWeights();
        calculateTurbulenceFields(turbModel, G(true), omega(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& omega0 = this->omega();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G =
        const_cast<FieldType&>
        (
            db().lookupObject<FieldType>(turbModel.GName())
        );

    FieldType& omega = const_cast<FieldType&>(internalField());

    forAll(*this, facei) {
        label celli = patch().faceCells()[facei];

        G[celli] = G0[celli];
        omega[celli] = omega0[celli];
    }

    fvPatchField<scalar>::updateCoeffs();
}


void omegaColebrookWallFunctionFvPatchScalarField::updateWeightedCoeffs
(
    const scalarField& weights
)
{
    if (updated()) {
        return;
    }
    Info<< "updateWeightedCoeffs" << nl << endl;

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
                                       (
                                           IOobject::groupName
                                           (
                                                   turbulenceModel::propertiesName,
                                                   internalField().group()
                                           )
                                       );

    setMaster();

    if (patch().index() == master_) {
        createAveragingWeights();
        calculateTurbulenceFields(turbModel, G(true), omega(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& omega0 = this->omega();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G =
        const_cast<FieldType&>
        (
            db().lookupObject<FieldType>(turbModel.GName())
        );

    FieldType& omega = const_cast<FieldType&>(internalField());

    scalarField& omegaf = *this;

    // only set the values if the weights are > tolerance
    forAll(weights, facei) {
        scalar w = weights[facei];

        if (w > tolerance_) {
            label celli = patch().faceCells()[facei];

            G[celli] = (1.0 - w)*G[celli] + w*G0[celli];
            omega[celli] = (1.0 - w)*omega[celli] + w*omega0[celli];
            omegaf[facei] = omega[celli];
        }
    }

    fvPatchField<scalar>::updateCoeffs();
}


void omegaColebrookWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix
)
{
    if (manipulatedMatrix()) {
        return;
    }

    matrix.setValues(patch().faceCells(), patchInternalField());

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void omegaColebrookWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix,
    const Field<scalar>& weights
)
{
    if (manipulatedMatrix()) {
        return;
    }

    DynamicList<label> constraintCells(weights.size());
    DynamicList<scalar> constraintomega(weights.size());
    const labelUList& faceCells = patch().faceCells();

    const DimensionedField<scalar, volMesh>& omega
        = internalField();

    label nConstrainedCells = 0;


    forAll(weights, facei) {
        // only set the values if the weights are > tolerance
        if (weights[facei] > tolerance_) {
            nConstrainedCells++;

            label celli = faceCells[facei];

            constraintCells.append(celli);
            constraintomega.append(omega[celli]);
        }
    }

    if (debug) {
        Pout<< "Patch: " << patch().name()
            << ": number of constrained cells = " << nConstrainedCells
            << " out of " << patch().size()
            << endl;
    }

    matrix.setValues
    (
        constraintCells,
        scalarField(constraintomega.xfer())
    );

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void omegaColebrookWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    writeLocalEntries(os);
    fixedValueFvPatchField<scalar>::write(os);
    os.writeKeyword("Betak") << Betak_ << token::END_STATEMENT << nl;
    os.writeKeyword("kr") << kr_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    omegaColebrookWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
