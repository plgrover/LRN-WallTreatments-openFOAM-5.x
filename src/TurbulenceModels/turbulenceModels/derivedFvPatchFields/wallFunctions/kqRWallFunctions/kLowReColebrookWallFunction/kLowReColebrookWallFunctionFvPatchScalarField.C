/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
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

#include "kLowReColebrookWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void kLowReColebrookWallFunctionFvPatchScalarField::checkType()
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


scalar kLowReColebrookWallFunctionFvPatchScalarField::yPlusLam
(
    const scalar kappa,
    const scalar E
)
{
    scalar ypl = 11.0;

    for (int i=0; i<10; i++) {
        ypl = log(max(E*ypl, 1))/kappa;
    }

    return ypl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kLowReColebrookWallFunctionFvPatchScalarField::kLowReColebrookWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
    :
    fixedValueFvPatchField<scalar>(p, iF),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8),
    Ceps2_(1.9),
    yPlusLam_(yPlusLam(kappa_, E_)),
    Betak_(0.09),
    kr_(0.002)
{
    checkType();
}


kLowReColebrookWallFunctionFvPatchScalarField::kLowReColebrookWallFunctionFvPatchScalarField
(
    const kLowReColebrookWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
    :
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_),
    Ceps2_(ptf.Ceps2_),
    yPlusLam_(ptf.yPlusLam_),
    Betak_(ptf.Betak_),
    kr_(ptf.kr_)
{
    checkType();
}


kLowReColebrookWallFunctionFvPatchScalarField::kLowReColebrookWallFunctionFvPatchScalarField
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
    Ceps2_(dict.lookupOrDefault<scalar>("Ceps2", 1.9)),
    yPlusLam_(yPlusLam(kappa_, E_)),
    Betak_(dict.lookupOrDefault<scalar>("Betak", 0.09)),
    kr_(dict.lookupOrDefault<scalar>("kr",0.002))
{
    checkType();
}


kLowReColebrookWallFunctionFvPatchScalarField::kLowReColebrookWallFunctionFvPatchScalarField
(
    const kLowReColebrookWallFunctionFvPatchScalarField& kwfpsf
)
    :
    fixedValueFvPatchField<scalar>(kwfpsf),
    Cmu_(kwfpsf.Cmu_),
    kappa_(kwfpsf.kappa_),
    E_(kwfpsf.E_),
    Ceps2_(kwfpsf.Ceps2_),
    yPlusLam_(kwfpsf.yPlusLam_),
    Betak_(kwfpsf.Betak_),
    kr_(kwfpsf.kr_)
{
    checkType();
}


kLowReColebrookWallFunctionFvPatchScalarField::kLowReColebrookWallFunctionFvPatchScalarField
(
    const kLowReColebrookWallFunctionFvPatchScalarField& kwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
    :
    fixedValueFvPatchField<scalar>(kwfpsf, iF),
    Cmu_(kwfpsf.Cmu_),
    kappa_(kwfpsf.kappa_),
    E_(kwfpsf.E_),
    Ceps2_(kwfpsf.Ceps2_),
    yPlusLam_(kwfpsf.yPlusLam_),
    Betak_(kwfpsf.Betak_),
    kr_(kwfpsf.kr_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void kLowReColebrookWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated()) {
        return;
    }

    const label patchi = patch().index();

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
                                    (
                                        IOobject::groupName
                                        (
                                            turbulenceModel::propertiesName,
                                            internalField().group()
                                        )
                                    );
    const scalarField& y = turbModel.y()[patchi];

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    // Get the kinematic and turbulent viscosities
    const tmp<volScalarField> tnu = turbModel.nu();
    const scalarField& nuw = tnu().boundaryField()[patchi];

    // Get the turbulent viscosity at the wall
    const tmp<volScalarField> tnut = turbModel.nut();
    const volScalarField& nut = tnut();
    const scalarField& nutw = nut.boundaryField()[patchi];

    // Get the velocity gradient at the wall
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magGradUw(mag(Uw.snGrad()));

    scalarField& kw = *this;

    // Set k wall values
    forAll(kw, facei) {
        label celli = patch().faceCells()[facei];
        scalar nueff = nutw[facei] + nuw[facei];
        scalar ustar = sqrt(nueff*magGradUw[facei]);

        if (ustar > ROOTVSMALL) {
            // Calculate kr+ - the non-dimensional roughness height
            // EQ (1) in Aupoix 2014
            scalar krPlus = kr_*ustar/nuw[facei];
            
            // Calculate the non-dimensional turbulent kinetic energy
            // Eq. 25 in Aupooix (2016)
            scalar kPlus = (1.0/sqrt(Betak_))
                *tanh(((log(krPlus/30.0)/log(10.)) + 1.0 
                - tanh(krPlus/125.))*tanh(krPlus/125.));

            kw[facei] = max(kPlus*(sqr(ustar)), 0.0);
        }
    }

    // Limit kw to avoid failure of the turbulence model due to division by kw
    kw = max(kw, SMALL);

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void kLowReColebrookWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    fixedValueFvPatchField<scalar>::evaluate(commsType);
}


void kLowReColebrookWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    os.writeKeyword("Cmu") << Cmu_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << kappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("E") << E_ << token::END_STATEMENT << nl;
    os.writeKeyword("Ceps2") << Ceps2_ << token::END_STATEMENT << nl;
    os.writeKeyword("Betak") << Betak_ << token::END_STATEMENT << nl;
    os.writeKeyword("kr") << kr_ << token::END_STATEMENT << nl;
    fixedValueFvPatchField<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    kLowReColebrookWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
