#ifndef staticPRiseAdjointPOutletFvPatchScalarField_H
#define staticPRiseAdjointPOutletFvPatchScalarField_H


#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class staticPRiseAdjointPOutletFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

public:

    //- Runtime type information
    TypeName("staticPRiseAdjointPOutlet");


    // Constructors

        //- Construct from patch and internal field
        staticPRiseAdjointPOutletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        staticPRiseAdjointPOutletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given staticPRiseAdjointPOutletFvPatchScalarField
        //  onto a new patch
        staticPRiseAdjointPOutletFvPatchScalarField
        (
            const staticPRiseAdjointPOutletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new staticPRiseAdjointPOutletFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        staticPRiseAdjointPOutletFvPatchScalarField
        (
            const staticPRiseAdjointPOutletFvPatchScalarField&,
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
                new staticPRiseAdjointPOutletFvPatchScalarField(*this, iF)
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
