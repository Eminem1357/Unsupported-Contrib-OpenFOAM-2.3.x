    const dictionary& smooth =
        mesh.solutionDict().subDict("smoothingControls");

    const int nSmoothingSteps =
        smooth.lookupOrDefault<int>("nSmoothingSteps", 0);
