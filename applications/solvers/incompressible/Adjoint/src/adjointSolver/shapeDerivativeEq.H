{
    Info<< nl << "Computing shape derivative" << endl;

    for 
    (
        int smoothingSteps=0; 
        smoothingSteps<=nSmoothingSteps; 
        smoothingSteps++
    )
    {

        fvVectorMatrix shapeDerivEqn
        (
            fvm::laplacian
            (
                dimensionedScalar
                (
                    "k",
                    dimensionSet(0, 0, 0, 0, 0),
                    0.001
                ),
                shapeDeriv
            )
        );

        solve(shapeDerivEqn);

    }
}

