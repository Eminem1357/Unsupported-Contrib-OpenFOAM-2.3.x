#ifndef staticPRiseAdjointUInletFvPatchVectorField_H
#define staticPRiseAdjointUInletFvPatchVectorField_H


#include "fixedValueFvPatchFields.H"

namespace Foam
{

class staticPRiseAdjointUInletFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

public:

    //- Runtime type information
    TypeName("staticPRiseAdjointUInlet");


    // Constructors

        //- Construct from patch and internal field
        staticPRiseAdjointUInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        staticPRiseAdjointUInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given staticPRiseAdjointUInletFvPatchVectorField
        //  onto a new patch
        staticPRiseAdjointUInletFvPatchVectorField
        (
            const staticPRiseAdjointUInletFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new staticPRiseAdjointUInletFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        staticPRiseAdjointUInletFvPatchVectorField
        (
            const staticPRiseAdjointUInletFvPatchVectorField&,
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
                new staticPRiseAdjointUInletFvPatchVectorField(*this, iF)
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
