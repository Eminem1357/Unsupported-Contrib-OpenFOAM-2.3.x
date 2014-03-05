
#include "staticPRiseAdjointUInletFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvPatchFieldMapper.H"

Foam::staticPRiseAdjointUInletFvPatchVectorField::
staticPRiseAdjointUInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


Foam::staticPRiseAdjointUInletFvPatchVectorField::
staticPRiseAdjointUInletFvPatchVectorField
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


Foam::staticPRiseAdjointUInletFvPatchVectorField::
staticPRiseAdjointUInletFvPatchVectorField
(
    const staticPRiseAdjointUInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::staticPRiseAdjointUInletFvPatchVectorField::
staticPRiseAdjointUInletFvPatchVectorField
(
    const staticPRiseAdjointUInletFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF)
{}


// Update the coefficients associated with the patch field
void Foam::staticPRiseAdjointUInletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<vector>& patchU =
        patch().lookupPatchField<volVectorField, vector>("U");

    scalarField patchUNormal(mag(patch().nf() & patchU));

    vectorField::operator=(-patchUNormal*patch().Sf()/(patch().magSf()));

    fixedValueFvPatchVectorField::updateCoeffs();
}


void 
Foam::staticPRiseAdjointUInletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry("value", os);
}


namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        staticPRiseAdjointUInletFvPatchVectorField
    );
}
