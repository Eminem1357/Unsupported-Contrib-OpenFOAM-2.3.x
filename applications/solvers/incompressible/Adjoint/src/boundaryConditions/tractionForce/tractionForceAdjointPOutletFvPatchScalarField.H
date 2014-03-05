#ifndef tractionForceAdjointPOutletFvPatchScalarField_H
#define tractionForceAdjointPOutletFvPatchScalarField_H


#include "fixedValueFvPatchFields.H"

namespace Foam
{

class tractionForceAdjointPOutletFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

public:

    //- Runtime type information
    TypeName("tractionForceAdjointPOutlet");


    // Constructors

        //- Construct from patch and internal field
        tractionForceAdjointPOutletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        tractionForceAdjointPOutletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given tractionForceAdjointPOutletFvPatchScalarField
        //  onto a new patch
        tractionForceAdjointPOutletFvPatchScalarField
        (
            const tractionForceAdjointPOutletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new tractionForceAdjointPOutletFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        tractionForceAdjointPOutletFvPatchScalarField
        (
            const tractionForceAdjointPOutletFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new tractionForceAdjointPOutletFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


} // End namespace Foam

#endif

