/*---------------------------------------------------------------------------*\
                pdfFoam: General Purpose PDF Solution Algorithm
                   for Reactive Flow Simulations in OpenFOAM

 Copyright (C) 2012 Michael Wild, Heng Xiao, Patrick Jenny,
                    Institute of Fluid Dynamics, ETH Zurich
-------------------------------------------------------------------------------
License
    This file is part of pdfFoam.

    pdfFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) version 3 of the same License.

    pdfFoam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with pdfFoam.  If not, see <http://www.gnu.org/licenses/>.

Application
    pdfSimpleFoam

Description
    Steady-state SIMPLE solver for laminar or turbulent RANS flow of
    compressible fluids.

Author
    Michael Wild

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "mcThermo.H"
#include "RASModel.H"
#include "simpleControl.H"
#include "turbulentFluidThermoModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    simpleControl simple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    bool prevCycleWasFV = false;
    volVectorField gradP = fvc::grad(p);

    scalar eqnResidual = 1, maxFVResidual = 0, maxPDFResidual = 0;
    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readThermoControls.H"

        if (FVCycle)
        {

            p.storePrevIter();
            rho.storePrevIter();

            // Pressure-velocity SIMPLE corrector
            {
                // TODO bound rho?
                rho.relax();
                Info<< "rho max/min : " << max(rho).value() << " " << min(rho).value() << endl;

                #include "UEqn.H"
                #include "pEqn.H"
            }
            turbulence->correct();
            prevCycleWasFV = true;
        }
        else
        {
            maxPDFResidual = 0.;

            if (prevCycleWasFV)
            {
                gradP = fvc::grad(p);
            }
            maxPDFResidual = thermo.evolve();
            prevCycleWasFV = false;
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
