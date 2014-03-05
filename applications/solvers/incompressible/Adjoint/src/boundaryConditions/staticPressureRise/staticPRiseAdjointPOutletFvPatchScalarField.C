#include "staticPRiseAdjointPOutletFvPatchScalarField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvPatchFieldMapper.H"
#include "fvc.H"

Foam::staticPRiseAdjointPOutletFvPatchScalarField::
staticPRiseAdjointPOutletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


Foam::staticPRiseAdjointPOutletFvPatchScalarField::
staticPRiseAdjointPOutletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF)
{
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));
}


Foam::staticPRiseAdjointPOutletFvPatchScalarField::
staticPRiseAdjointPOutletFvPatchScalarField
(
    const staticPRiseAdjointPOutletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::staticPRiseAdjointPOutletFvPatchScalarField::
staticPRiseAdjointPOutletFvPatchScalarField
(
    const staticPRiseAdjointPOutletFvPatchScalarField& pivpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(pivpvf, iF)
{}


// Update the coefficients associated with the patch field
void Foam::staticPRiseAdjointPOutletFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<vector>& patchU =
        patch().lookupPatchField<volVectorField, scalar>("U");

    scalarField patchUNormal(mag(patch().nf() & patchU));

    const fvPatchField<scalar>& patchP =
        patch().lookupPatchField<volScalarField, scalar>("p");

    const fvPatchField<vector>& patchAdjointU =
        patch().lookupPatchField<volVectorField, scalar>("adjointU");

    scalarField::operator=
    (
        (patchU & patchAdjointU) + patchP
    );

    fixedValueFvPatchScalarField::updateCoeffs();
}


void 
Foam::staticPRiseAdjointPOutletFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry("value", os);
}

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        staticPRiseAdjointPOutletFvPatchScalarField
    );
}
