#ifndef dissipatedPowerAdjointUInletFvPatchVectorField_H
#define dissipatedPowerAdjointUInletFvPatchVectorField_H


#include "fixedValueFvPatchFields.H"

namespace Foam
{

class dissipatedPowerAdjointUInletFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

public:

    //- Runtime type information
    TypeName("dissipatedPowerAdjointUInlet");


    // Constructors

        //- Construct from patch and internal field
        dissipatedPowerAdjointUInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        dissipatedPowerAdjointUInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given 
        //  dissipatedPowerAdjointUInletFvPatchVectorField
        //  onto a new patch
        dissipatedPowerAdjointUInletFvPatchVectorField
        (
            const dissipatedPowerAdjointUInletFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new dissipatedPowerAdjointUInletFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        dissipatedPowerAdjointUInletFvPatchVectorField
        (
            const dissipatedPowerAdjointUInletFvPatchVectorField&,
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
                new dissipatedPowerAdjointUInletFvPatchVectorField(*this, iF)
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
