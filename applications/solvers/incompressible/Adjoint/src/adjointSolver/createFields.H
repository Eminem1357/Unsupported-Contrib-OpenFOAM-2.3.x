    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar nu
    (
        transportProperties.lookup("nu")
    );

    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    surfaceScalarField phiPrime
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(U) & mesh.Sf()
    );

    adjustPhi(phiPrime, U, p);

    Info << "Reading field adjointU\n" << endl;
    volVectorField adjointU
    (
        IOobject
        (
            "adjointU",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field adjointP\n" << endl;
    volScalarField adjointP
    (
        IOobject
        (
            "adjointP",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    surfaceScalarField adjointPhi
    (
        IOobject
        (
            "adjointPhi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(adjointU) & mesh.Sf()
    );

    label adjointPRefCell = 0;
    scalar adjointPRefValue = 0.0;
    setRefCell
    (
        adjointP, 
        mesh.solutionDict().subDict("SIMPLE"), 
        adjointPRefCell, 
        adjointPRefValue
    );

    Info<< "Reading field nut\n" << endl;
    volScalarField nut
    (
        IOobject
        (
            "nut",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField nuEff
    (
        IOobject
        (
            "nuEff",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        nu + nut
    );

    surfaceVectorField n = mesh.Sf()/mag(mesh.Sf());

    surfaceVectorField faceCenters
    (
        IOobject
        (
            "faceCenters",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh.Cf()*scalar(1)
    );

    volVectorField cellCenters
    (
        IOobject
        (
            "cellCenters",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh.C()*scalar(1)
    );

    surfaceVectorField faceNormals
    (
        IOobject
        (
            "faceNormals",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        n*scalar(1)
    );

    surfaceVectorField s
    (
        IOobject
        (
            "s",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        (
	    ((linearInterpolate(U) & n) * (linearInterpolate(adjointU))) 
          +
            (linearInterpolate(adjointP) * n)
          +
            (
                 (linearInterpolate(nuEff*fvc::grad(adjointU)) & n)
               +
                 (n & linearInterpolate(nuEff*T(fvc::grad(adjointU))))
	    )
        )
    );


    surfaceVectorField sShapeDerivative
    (
        IOobject
        (
            "sShapeDerivative",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        -(s & (linearInterpolate(fvc::grad(U))))*mag(mesh.Sf())
    );

    surfaceVectorField sNormalShapeDerivative
    (
        IOobject
        (
            "sNormalShapeDerivative",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        -((s & (linearInterpolate(fvc::grad(U)))) & n)*mesh.Sf()
    );

    dimensionedVector initVec
    (
        "initVec",
        dimensionSet(0, 4, -3, 0, 0, 0, 0),
        vector(0, 0, 0)
    );

    volVectorField shapeDerivative
    (
        IOobject
        (
            "shapeDerivative",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        initVec
    );

    volVectorField normalShapeDerivative
    (
        IOobject
        (
            "normalShapeDerivative",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        initVec
    );


    volVectorField shapeDeriv
    (
        IOobject
        (
            "shapeDeriv",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
