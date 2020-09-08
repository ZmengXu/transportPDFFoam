    // Solve the Momentum equation

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
    );

    fvVectorMatrix& UEqn = tUEqn.ref();
    UEqn.relax();

    solve(UEqn == -fvc::grad(p));

    maxFVResidual = max(eqnResidual, maxFVResidual);
