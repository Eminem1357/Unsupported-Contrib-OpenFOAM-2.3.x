/* 
   Steady-state adjoint SIMPLE solver for turbulent flows.
   The adjoint formulation is based on the discrete adjoint formal
   relation that preserves time symmetry in continuous limit.
   This adjoint solver is stable in energy norm.

*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "RASModel.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readSmoothingControls.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    simpleControl simple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // adjoint SIMPLE predictor-corrector
        {
            #include "adjointUEqn.H"
            #include "adjointPEqn.H"
        }

        #include "shapeSensitivity.H"
        #include "shapeDerivative.H"

        if (runTime.write())
            #include "shapeDerivativeEq.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
