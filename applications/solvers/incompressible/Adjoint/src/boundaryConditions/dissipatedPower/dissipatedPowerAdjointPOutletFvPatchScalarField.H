#ifndef dissipatedPowerAdjointPOutletFvPatchScalarField_H
#define dissipatedPowerAdjointPOutletFvPatchScalarField_H


#include "fixedValueFvPatchFields.H"

namespace Foam
{

class dissipatedPowerAdjointPOutletFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{

public:

    //- Runtime type information
    TypeName("dissipatedPowerAdjointPOutlet");


    // Constructors

        //- Construct from patch and internal field
        dissipatedPowerAdjointPOutletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        dissipatedPowerAdjointPOutletFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given 
        //  dissipatedPowerAdjointPOutletFvPatchScalarField
        //  onto a new patch
        dissipatedPowerAdjointPOutletFvPatchScalarField
        (
            const dissipatedPowerAdjointPOutletFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new dissipatedPowerAdjointPOutletFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        dissipatedPowerAdjointPOutletFvPatchScalarField
        (
            const dissipatedPowerAdjointPOutletFvPatchScalarField&,
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
                new dissipatedPowerAdjointPOutletFvPatchScalarField(*this, iF)
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
