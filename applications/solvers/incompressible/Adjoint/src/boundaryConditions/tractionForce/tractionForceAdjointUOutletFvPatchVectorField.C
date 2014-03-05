#include "tractionForceAdjointUOutletFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvPatchFieldMapper.H"

Foam::tractionForceAdjointUOutletFvPatchVectorField::
tractionForceAdjointUOutletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


Foam::tractionForceAdjointUOutletFvPatchVectorField::
tractionForceAdjointUOutletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF)
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::tractionForceAdjointUOutletFvPatchVectorField::
tractionForceAdjointUOutletFvPatchVectorField
(
    const tractionForceAdjointUOutletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::tractionForceAdjointUOutletFvPatchVectorField::
tractionForceAdjointUOutletFvPatchVectorField
(
    const tractionForceAdjointUOutletFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF)
{}


// Update the coefficients associated with the patch field
void Foam::tractionForceAdjointUOutletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<vector>& patchU =
        patch().lookupPatchField<volVectorField, vector>("U");
    const fvPatchField<scalar>& patchAdjointP =
        patch().lookupPatchField<volScalarField, scalar>("adjointP");

    scalarField patchUNormal(mag(patch().nf() & patchU));

    vectorField::operator=
        ((-patchAdjointP/patchUNormal)*patch().nf());

    fixedValueFvPatchVectorField::updateCoeffs();
}


void 
Foam::tractionForceAdjointUOutletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);
}


namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        tractionForceAdjointUOutletFvPatchVectorField
    );
}
