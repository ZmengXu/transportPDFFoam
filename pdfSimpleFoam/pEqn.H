volScalarField rUA = 1.0/UEqn.A();
U = rUA*UEqn.H();
tUEqn.clear();

phi = fvc::interpolate(rho)*(fvc::interpolate(U) & mesh.Sf());

while (simple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rho*rUA, p) == fvc::div(phi)
    );

    pEqn.setReference(pRefCell, pRefValue);

    pEqn.solve();

    if (simple.finalNonOrthogonalIter())
    {
        phi -= pEqn.flux();
    }
}


#include "incompressible/continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

volVectorField gradP = fvc::grad(p);
U -= rUA*gradP;
U.correctBoundaryConditions();
