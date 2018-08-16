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

Class
    Foam::kLowReNikuradseWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides a turbulence kinetic energy wall function 
    condition for low-Reynolds number turbulent flow cases over rough walls 
    using the Nikuradse based correction.

    The rough wall correction is described in:
    \verbatim
	    Aupoix, B. (2015). Roughness corrections for the k–ω shear stress 
        transport model: Status and proposals. Journal of Fluids Engineering, 
        137(2), 021202.
    \endverbatim


Usage
    \table
        Property     | Description             | Required    | Default value
        Cmu          | model coefficient       | no          | 0.09
        kappa        | Von Karman constant     | no          | 0.41
        E            | model coefficient       | no          | 9.8
        Ceps2        | model coefficient       | no          | 1.9
	    kr	         | sand grain height       | no          | 0.002
	    Betak        | model coefficient       |             | 0.09
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            kLowReNikuradseWallFunction;
    }
    \endverbatim

See also
    Foam::kLowReWallFunctionFvPatchScalarField

SourceFiles
    kLowReNikuradseWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef kLowReNikuradseWallFunctionFvPatchScalarField_H
#define kLowReNikuradseWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class kLowReNikuradseWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class kLowReNikuradseWallFunctionFvPatchScalarField
    :
public fixedValueFvPatchField<scalar>
{
protected:

    // Protected data

    //- Cmu coefficient
    scalar Cmu_;

    //- Von Karman constant
    scalar kappa_;

    //- E coefficient
    scalar E_;

    //- Ceps2 coefficient
    scalar Ceps2_;

    //- Y+ at the edge of the laminar sublayer
    scalar yPlusLam_;

    //- kr Nikuradse roughness coefficient
    scalar kr_;

    //- Betak_ constant of the k-ω model
    scalar Betak_;

    // Protected Member Functions

    //- Check the type of the patch
    virtual void checkType();

    //- Calculate the Y+ at the edge of the laminar sublayer
    scalar yPlusLam(const scalar kappa, const scalar E);


public:

    //- Runtime type information
    TypeName("kLowReNikuradseWallFunction");


    // Constructors

    //- Construct from patch and internal field
    kLowReNikuradseWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct from patch, internal field and dictionary
    kLowReNikuradseWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    //- Construct by mapping given kLowReNikuradseWallFunctionFvPatchScalarField
    //  onto a new patch
    kLowReNikuradseWallFunctionFvPatchScalarField
    (
        const kLowReNikuradseWallFunctionFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    //- Construct as copy
    kLowReNikuradseWallFunctionFvPatchScalarField
    (
        const kLowReNikuradseWallFunctionFvPatchScalarField&
    );

    //- Construct and return a clone
    virtual tmp<fvPatchScalarField> clone() const {
        return tmp<fvPatchScalarField>
               (
                   new kLowReNikuradseWallFunctionFvPatchScalarField(*this)
               );
    }

    //- Construct as copy setting internal field reference
    kLowReNikuradseWallFunctionFvPatchScalarField
    (
        const kLowReNikuradseWallFunctionFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const {
        return tmp<fvPatchScalarField>
               (
                   new kLowReNikuradseWallFunctionFvPatchScalarField(*this, iF)
               );
    }

    // Member functions

    // Evaluation functions

    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();

    //- Evaluate the patchField
    virtual void evaluate(const Pstream::commsTypes);


    // I-O

    //- Write
    virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
