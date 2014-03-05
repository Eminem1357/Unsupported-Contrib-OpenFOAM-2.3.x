#ifndef tractionForceAdjointUOutletFvPatchVectorField_H
#define tractionForceAdjointUOutletFvPatchVectorField_H


#include "fixedValueFvPatchFields.H"

namespace Foam
{

class tractionForceAdjointUOutletFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

public:

    //- Runtime type information
    TypeName("tractionForceAdjointUOutlet");


    // Constructors

        //- Construct from patch and internal field
        tractionForceAdjointUOutletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        tractionForceAdjointUOutletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given tractionForceAdjointUOutletFvPatchVectorField
        //  onto a new patch
        tractionForceAdjointUOutletFvPatchVectorField
        (
            const tractionForceAdjointUOutletFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new tractionForceAdjointUOutletFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        tractionForceAdjointUOutletFvPatchVectorField
        (
            const tractionForceAdjointUOutletFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new tractionForceAdjointUOutletFvPatchVectorField(*this, iF)
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
