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

Class
    Foam::omegaNikuradseWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides a turbulence specific dissipation, omega,
    wall function condition for low-Reynolds number turbulent flow cases over 
    rough walls using the Nikuradse based correction.  
     
    The rough wall correction is described in:
    \verbatim
	    Aupoix, B. (2015). Roughness corrections for the k–ω shear stress 
        transport model: Status and proposals. Journal of Fluids Engineering, 
        137(2), 021202.
    \endverbatim


Usage
    \table
        Property     | Description             | Required    | Default value
        Cmu          | Model coefficient       | no          | 0.09
        kappa        | von Karman constant     | no          | 0.41
        E            | Model coefficient       | no          | 9.8
        beta1        | Model coefficient       | no          | 0.075
        blended      | Blending switch         | no          | false
	    kr	         | sand grain height       | no          | 0.002
	    Betak        | model coefficient       |             | 0.09
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            omegaNikuradseWallFunction;
    }
    \endverbatim

See also
    Foam::fixedInternalValueFvPatchField
    Foam::omegaWallFunctionFvPatchScalarField
    Foam::epsilonWallFunctionFvPatchScalarField

SourceFiles
    omegaNikuradseWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef omegaNikuradseWallFunctionFvPatchScalarField_H
#define omegaNikuradseWallFunctionFvPatchScalarField_H

#include "fixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class turbulenceModel;

/*---------------------------------------------------------------------------*\
             Class omegaNikuradseWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class omegaNikuradseWallFunctionFvPatchScalarField
    :
public fixedValueFvPatchField<scalar>
{
protected:

    // Protected data

    //- Tolerance used in weighted calculations
    static scalar tolerance_;

    //- Cmu coefficient
    scalar Cmu_;

    //- Von Karman constant
    scalar kappa_;

    //- E coefficient
    scalar E_;

    //- beta1 coefficient
    scalar beta1_;

    //- beta1 coefficient
    Switch blended_;

    //- y+ at the edge of the laminar sublayer
    scalar yPlusLam_;

    //- kr Nikuradse roughness coefficient
    scalar kr_;

    //- Betak_constant of the k-ω model
    scalar Betak_;

    //- Local copy of turbulence G field
    scalarField G_;

    //- Local copy of turbulence omega field
    scalarField omega_;

    //- Initialised flag
    bool initialised_;

    //- Master patch ID
    label master_;

    //- List of averaging corner weights
    List<List<scalar>> cornerWeights_;


    // Protected Member Functions

    //- Check the type of the patch
    virtual void checkType();

    //- Write local wall function variables
    virtual void writeLocalEntries(Ostream&) const;

    //- Set the master patch - master is responsible for updating all
    //  wall function patches
    virtual void setMaster();

    //- Create the averaging weights for cells which are bounded by
    //  multiple wall function faces
    virtual void createAveragingWeights();

    //- Helper function to return non-const access to an omega patch
    virtual omegaNikuradseWallFunctionFvPatchScalarField& omegaPatch
    (
        const label patchi
    );

    //- Main driver to calculate the turbulence fields
    virtual void calculateTurbulenceFields
    (
        const turbulenceModel& turbulence,
        scalarField& G0,
        scalarField& omega0
    );

    //- Calculate the omega and G
    virtual void calculate
    (
        const turbulenceModel& turbulence,
        const List<scalar>& cornerWeights,
        const fvPatch& patch,
        scalarField& G,
        scalarField& omega
    );

    //- Return non-const access to the master patch ID
    virtual label& master() {
        return master_;
    }


public:

    //- Runtime type information
    TypeName("omegaNikuradseWallFunction");


    // Constructors

    //- Construct from patch and internal field
    omegaNikuradseWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct from patch, internal field and dictionary
    omegaNikuradseWallFunctionFvPatchScalarField
    (
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const dictionary&
    );

    //- Construct by mapping given
    // omegaNikuradseWallFunctionFvPatchScalarField
    //  onto a new patch
    omegaNikuradseWallFunctionFvPatchScalarField
    (
        const omegaNikuradseWallFunctionFvPatchScalarField&,
        const fvPatch&,
        const DimensionedField<scalar, volMesh>&,
        const fvPatchFieldMapper&
    );

    //- Construct as copy
    omegaNikuradseWallFunctionFvPatchScalarField
    (
        const omegaNikuradseWallFunctionFvPatchScalarField&
    );

    //- Construct and return a clone
    virtual tmp<fvPatchScalarField> clone() const {
        return tmp<fvPatchScalarField>
               (
                   new omegaNikuradseWallFunctionFvPatchScalarField(*this)
               );
    }

    //- Construct as copy setting internal field reference
    omegaNikuradseWallFunctionFvPatchScalarField
    (
        const omegaNikuradseWallFunctionFvPatchScalarField&,
        const DimensionedField<scalar, volMesh>&
    );

    //- Construct and return a clone setting internal field reference
    virtual tmp<fvPatchScalarField> clone
    (
        const DimensionedField<scalar, volMesh>& iF
    ) const {
        return tmp<fvPatchScalarField>
               (
                   new omegaNikuradseWallFunctionFvPatchScalarField(*this, iF)
               );
    }


    // Member functions

    // Access

    //- Return non-const access to the master's G field
    scalarField& G(bool init = false);

    //- Return non-const access to the master's omega field
    scalarField& omega(bool init = false);


    // Evaluation functions

    //- Update the coefficients associated with the patch field
    virtual void updateCoeffs();

    //- Update the coefficients associated with the patch field
    virtual void updateWeightedCoeffs(const scalarField& weights);

    //- Manipulate matrix
    virtual void manipulateMatrix(fvMatrix<scalar>& matrix);

    //- Manipulate matrix with given weights
    virtual void manipulateMatrix
    (
        fvMatrix<scalar>& matrix,
        const scalarField& weights
    );


    // I-O

    //- Write
    virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
