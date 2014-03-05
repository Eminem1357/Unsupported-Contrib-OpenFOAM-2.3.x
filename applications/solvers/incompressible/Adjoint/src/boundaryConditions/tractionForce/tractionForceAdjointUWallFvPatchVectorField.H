#ifndef tractionForceAdjointUWallFvPatchVectorField_H
#define tractionForceAdjointUWallFvPatchVectorField_H


#include "fixedValueFvPatchFields.H"

namespace Foam
{

class tractionForceAdjointUWallFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{

public:

    //- Runtime type information
    TypeName("tractionForceAdjointUWall");


    // Constructors

        //- Construct from patch and internal field
        tractionForceAdjointUWallFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        tractionForceAdjointUWallFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given 
        //  tractionForceAdjointUWallFvPatchVectorField
        //  onto a new patch
        tractionForceAdjointUWallFvPatchVectorField
        (
            const tractionForceAdjointUWallFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new tractionForceAdjointUWallFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        tractionForceAdjointUWallFvPatchVectorField
        (
            const tractionForceAdjointUWallFvPatchVectorField&,
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
                new tractionForceAdjointUWallFvPatchVectorField(*this, iF)
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


