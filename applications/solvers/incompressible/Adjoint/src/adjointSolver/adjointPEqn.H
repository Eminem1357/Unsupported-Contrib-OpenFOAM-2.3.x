// Adjoint pressure equation and corrector
{
    adjointP.boundaryField().updateCoeffs();

    volScalarField rAadjointU(1.0/adjointUEqn().A());
    adjointU = rAadjointU*adjointUEqn().H();
    adjointUEqn.clear();

    adjointPhi = fvc::interpolate(adjointU, "interpolate(HbyA)") & mesh.Sf();
    adjustPhi(adjointPhi, adjointU, adjointP);

    // Non-orthogonal pressure corrector loop
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix adjointPEqn
        (
            fvm::laplacian(rAadjointU, adjointP) == fvc::div(adjointPhi)
        );

        adjointPEqn.setReference(adjointPRefCell, adjointPRefValue);

        adjointPEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            adjointPhi -= adjointPEqn.flux();
        }
    }

    // Check for adjoint continuity errors
    {
        volScalarField contErr = fvc::div(adjointPhi);

        scalar sumLocalContErr = runTime.deltaT().value()*
            mag(contErr)().weightedAverage(mesh.V()).value();

        scalar globalContErr = runTime.deltaT().value()*
            contErr.weightedAverage(mesh.V()).value();
        cumulativeContErr += globalContErr;

        Info<< "time step adjoint continuity errors : sum local = " 
            << sumLocalContErr
            << ", global = " << globalContErr
            << ", cumulative = " << cumulativeContErr
            << endl;
    }


    // Explicitly relax pressure for momentum corrector
    adjointP.relax();

    // Adjoint momentum corrector
    adjointU -= rAadjointU*fvc::grad(adjointP);
    adjointU.correctBoundaryConditions();
}
