
#include "dissipatedPowerAdjointPOutletFvPatchScalarField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "surfaceFields.H"
#include "fvPatchFieldMapper.H"
#include "fvc.H"

Foam::dissipatedPowerAdjointPOutletFvPatchScalarField::
dissipatedPowerAdjointPOutletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF)
{}


Foam::dissipatedPowerAdjointPOutletFvPatchScalarField::
dissipatedPowerAdjointPOutletFvPatchScalarField
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


Foam::dissipatedPowerAdjointPOutletFvPatchScalarField::
dissipatedPowerAdjointPOutletFvPatchScalarField
(
    const dissipatedPowerAdjointPOutletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper)
{}


Foam::dissipatedPowerAdjointPOutletFvPatchScalarField::
dissipatedPowerAdjointPOutletFvPatchScalarField
(
    const dissipatedPowerAdjointPOutletFvPatchScalarField& pivpvf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(pivpvf, iF)
{}


// Update the coefficients associated with the patch field
void Foam::dissipatedPowerAdjointPOutletFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<vector>& patchU =
        patch().lookupPatchField<volVectorField, scalar>("U");

    scalarField patchUNormal(mag(patch().nf() & patchU));

    const fvPatchField<vector>& patchAdjointU =
        patch().lookupPatchField<volVectorField, scalar>("adjointU");

    scalarField patchAdjointUNormal(mag(patch().nf() & patchAdjointU));

    scalarField::operator=
    (
        (patchU & patchAdjointU) 
      - 0.5*magSqr(patchU) - magSqr(patchUNormal)
    );

    fixedValueFvPatchScalarField::updateCoeffs();
}


void 
Foam::dissipatedPowerAdjointPOutletFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry("value", os);
}

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        dissipatedPowerAdjointPOutletFvPatchScalarField
    );
}
