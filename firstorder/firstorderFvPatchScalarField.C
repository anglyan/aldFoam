/*---------------------------------------------------------------------------

Copyright 2019-2023 Argonne UChicago LLC

This file is part of aldFoam.

    aldFoam is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    aldFoam is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with aldFoam.  If not, see <https://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "firstorderFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::firstorderFvPatchScalarField::t() const
{
    return db().time().userTimeValue();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::firstorderFvPatchScalarField::
firstorderFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    fieldData_(p.size(), Zero),
    betaName_("beta"),
    diffName_("D1"),
    vthName_("vth")
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = 0.0;
}


Foam::firstorderFvPatchScalarField::
firstorderFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    fieldData_(p.size(), Zero),
    betaName_(dict.lookupOrDefault<word>("betaField", "beta")),
    diffName_(dict.lookupOrDefault<word>("diffCoeff", "D1")),
    vthName_(dict.lookupOrDefault<word>("vth", "vth"))
{
    refGrad() = Zero;
    valueFraction() = 0.0;

    refValue() = Zero;
    fvPatchScalarField::operator=(refValue());

}


Foam::firstorderFvPatchScalarField::
firstorderFvPatchScalarField
(
    const firstorderFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    fieldData_(mapper(ptf.fieldData_)),
    betaName_(ptf.betaName_),
    diffName_(ptf.diffName_),
    vthName_(ptf.vthName_)
{}


Foam::firstorderFvPatchScalarField::
firstorderFvPatchScalarField
(
    const firstorderFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    fieldData_(ptf.fieldData_),
    betaName_(ptf.betaName_),
    diffName_(ptf.diffName_),
    vthName_(ptf.vthName_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::firstorderFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
    m(fieldData_, fieldData_);
}


void Foam::firstorderFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const firstorderFvPatchScalarField& tiptf =
        refCast<const firstorderFvPatchScalarField>(ptf);

    fieldData_.rmap(tiptf.fieldData_, addr);
}


void Foam::firstorderFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    mixedFvPatchScalarField::reset(ptf);

    const firstorderFvPatchScalarField& tiptf =
        refCast<const firstorderFvPatchScalarField>(ptf);

    fieldData_.reset(tiptf.fieldData_);
}


void Foam::firstorderFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarField& beta =
        patch().template lookupPatchField<volScalarField, scalar>
            ( betaName_ );

    const dictionary& transportProperties=
        this->db().template lookupObject<IOdictionary>
            ("transportProperties");

    dimensionedScalar DN(transportProperties.lookup(diffName_));
    dimensionedScalar vth(transportProperties.lookup(vthName_));

    fieldData_ = 0.25*beta*vth.value()/DN.value();

    refGrad() = Zero;
    refValue() = Zero;
    valueFraction() = fieldData_/(fieldData_ + patch().deltaCoeffs());

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::firstorderFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        firstorderFvPatchScalarField
    );
}

// ************************************************************************* //
