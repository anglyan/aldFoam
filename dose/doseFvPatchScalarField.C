/*---------------------------------------------------------------------------

Copyright 2019-2025 Argonne UChicago LLC

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

#include "doseFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::doseFvPatchScalarField::t() const
{
    return db().time().value();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::doseFvPatchScalarField::
doseFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    doseValue("doseValue", dict, p.size()),
    doseStart(readScalar(dict.lookup("doseStart"))),
    doseEnd(readScalar(dict.lookup("doseEnd"))),
    diffName_(dict.lookupOrDefault<word>("diffCoeff","D1"))
{
    refGrad() = Zero;
    valueFraction() = 0.0;
    refValue() = doseValue;
    fvPatchScalarField::operator=(refValue());

}


Foam::doseFvPatchScalarField::
doseFvPatchScalarField
(
    const doseFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    doseValue(mapper(ptf.doseValue)),
    doseStart(ptf.doseStart),
    doseEnd(ptf.doseEnd),
    diffName_(ptf.diffName_)
{}


Foam::doseFvPatchScalarField::
doseFvPatchScalarField
(
    const doseFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    doseValue(ptf.doseValue),
    doseStart(ptf.doseStart),
    doseEnd(ptf.doseEnd),
    diffName_(ptf.diffName_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::doseFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fieldMapper& mapper
)
{
    mixedFvPatchScalarField::map(ptf, mapper);

    const doseFvPatchScalarField& tiptf =
        refCast<const doseFvPatchScalarField>(ptf);

    mapper(doseValue, tiptf.doseValue);
}


void Foam::doseFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    mixedFvPatchScalarField::reset(ptf);

    const doseFvPatchScalarField& tiptf =
        refCast<const doseFvPatchScalarField>(ptf);

    doseValue.reset(tiptf.doseValue);
}


void Foam::doseFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    const scalar currtime = t();

    if (currtime >= doseStart && currtime < doseEnd){
        refValue() = doseValue;
        refGrad() = 0.0;
        valueFraction() = 1.0;
    }
    else {
        const dictionary & transportProperties =
            this->db().lookupObject<IOdictionary>
                (
                    "transportProperties"
                );
        dimensionedScalar DT(transportProperties.lookup(diffName_));
        scalar diffCoeff(DT.value());
        const vectorField & Up =
            patch().lookupPatchField<volVectorField,vector>("U");
        scalarField Un =  -Up & patch().nf()/diffCoeff;
        
        valueFraction() = Un/(Un + patch().deltaCoeffs());
        refGrad() = 0.0;
        refValue() = 0.0;
//        valueFraction() = 0.0;
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::doseFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "doseValue", doseValue);
    os.writeKeyword("doseStart") << doseStart << token::END_STATEMENT <<nl;
    os.writeKeyword("doseEnd") << doseEnd << token::END_STATEMENT <<nl;
    writeEntry(os, "value", *this);
}



// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        doseFvPatchScalarField
    );
}

// ************************************************************************* //
