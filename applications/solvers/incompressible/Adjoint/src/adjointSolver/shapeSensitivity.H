    s = (
            ((linearInterpolate(U*adjointU) & n))
          +
            (linearInterpolate(adjointP) * n)
          +
            (
                 (linearInterpolate(nuEff*fvc::grad(adjointU)) & n)
               +
                 (n & linearInterpolate(nuEff*T(fvc::grad(adjointU))))
	    )
        );

    sShapeDerivative = -(s&(linearInterpolate(T(fvc::grad(U))))*mag(mesh.Sf()));

    sNormalShapeDerivative = (sShapeDerivative & n)*n;

    forAll(shapeDerivative.boundaryField(), patchi)
    {
        vectorField&   sdP = shapeDerivative.boundaryField()[patchi];
        vectorField&  ssdP = sShapeDerivative.boundaryField()[patchi];
        vectorField&   sn  = normalShapeDerivative.boundaryField()[patchi];
        vectorField&  ssn  = sNormalShapeDerivative.boundaryField()[patchi];
        forAll(sdP,facei)
        {
            sdP[facei] = ssdP[facei];
             sn[facei] = ssn[facei];
        }
#if 0
        const scalar alpha = 0.1396263;//8 degrees //0.03410651;// 1 degree
        vectorField&  rotS = sRotShapeDerivative.boundaryField()[patchi];
        vectorField& coord = faceCenters.boundaryField()[patchi];
        forAll(rotS,facei)
        {
	  rotS[facei][0] = 
              ssdP[facei][0]*(coord[facei][0] - 3.165)*(Foam::cos(alpha)-1.0) 
            - ssdP[facei][2]*(coord[facei][2] - 0.803)*Foam::sin(alpha);
          rotS[facei][2] = 
              ssdP[facei][0]*(coord[facei][0] - 3.165)*Foam::sin(alpha) 
            - ssdP[facei][2]*(coord[facei][2] - 0.803)*(Foam::cos(alpha)-1.0);
          rotS[facei][1] = ssdP[facei][1]; 
        }

#endif
    }
    


