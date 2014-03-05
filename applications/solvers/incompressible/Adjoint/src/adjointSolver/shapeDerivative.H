    const fvPatchList& patches = mesh.boundary();

    forAll(patches, patchI)
    {
        const fvPatch& currPatch = patches[patchI];
        fvsPatchVectorField& pS = sShapeDerivative.boundaryField()[patchI];
        forAll(currPatch, faceI)
        {
            label i = currPatch.faceCells()[faceI];
            shapeDeriv[i] = pS[faceI];
        }
    }


