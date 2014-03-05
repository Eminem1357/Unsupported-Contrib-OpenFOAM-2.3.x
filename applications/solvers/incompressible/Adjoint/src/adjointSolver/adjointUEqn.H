    // Adjoint momentum predictor

    // Note: adjoint momentum equation consistent with the kinetic
    // energy conservation (time symmetry enforced by the selection
    // of the linear form of the momentum operator)

    tmp<fvVectorMatrix> adjointUEqn
    (
        fvm::div(-phiPrime,adjointU)
      - fvm::laplacian(nuEff,adjointU)
      - fvc::div(-nuEff*dev(T(-fvc::grad(adjointU))))
    );

    adjointUEqn().relax();

    solve(adjointUEqn() == -(-fvc::grad(adjointP)));

