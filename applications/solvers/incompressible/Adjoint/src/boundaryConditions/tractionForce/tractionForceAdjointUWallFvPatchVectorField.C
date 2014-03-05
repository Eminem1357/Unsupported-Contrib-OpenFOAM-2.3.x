
#include "tractionForceAdjointUWallFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvPatchFieldMapper.H"

Foam::tractionForceAdjointUWallFvPatchVectorField::
tractionForceAdjointUWallFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


Foam::tractionForceAdjointUWallFvPatchVectorField::
tractionForceAdjointUWallFvPatchVectorField
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


Foam::tractionForceAdjointUWallFvPatchVectorField::
tractionForceAdjointUWallFvPatchVectorField
(
    const tractionForceAdjointUWallFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::tractionForceAdjointUWallFvPatchVectorField::
tractionForceAdjointUWallFvPatchVectorField
(
    const tractionForceAdjointUWallFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF)
{}


// Update the coefficients associated with the patch field
void Foam::tractionForceAdjointUWallFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    vectorField::operator=(-patch().Sf()/(patch().magSf()));

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::tractionForceAdjointUWallFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);
}


namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        tractionForceAdjointUWallFvPatchVectorField
    );
}
