/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
PISOFoam - Pressure Implicit with Splitting of Operators
Transient Solver for Incompressible Flow
Turbulence modelling is generic (laminar, RAS, LES)
Algorithm main steps:
    1) Set the boundary condition
    2) Solve the discretized momenum equation to compute an intermediate velocity field
    3) Compute the mass fluxes at the cell faces
    4) Solve the pressure equation
    5) Correct the mass fluxes at the cell faces
    6) Correct the velocities on the basis of the new pressure field
    7) Update the boundary conditions
    8) Step 3 - 7 are repeated for the prescribed number of correctors
    9) Increase the time step and repeat from step 2

    KALMAN FILTER - PARALLEL VERSION - REAL TIME MEASUREMENT
    

(Uncleaned solver that is a copy of this one is pisoFoamTestFlag!!!!)


\*----------------------------------------------------------------------------------------------------------------------------*/

// - Standard CFD Library (always included)
// - Class definition of the single-phase transport model based on viscosity
// - Abstract base class for turbulence models (RAS, LES, laminar)
// - Library for handling traditional matrices
// - Library for reading from external files (Input)

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "scalarMatrices.H"
#include "IFstream.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    // - Check root path and case path
    // - Creates object runTime of class Time to control time during the simulation
    // - Creates the object mesh of class fvMesh (which contains all information related to the mesh)
    // - Declare the cumulative, local and global continuity errors (all initialized at zero)
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"
#   include "createFields.H"
#   include "initContinuityErrs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
#   include "getProperties.H"
#   include "initialize.H"
    

    
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        

        // - Reads from the PISO subdictionary in the fvSolution file
        // nCorr - Integer containing the number of correction loops
        // nNonOrthCorr - Integer containing the number of non orthogonal correctors (default 0)
        // momentumPredictor - Boolean that indicates whether the velocity equation is solved (default TRUE)
        // transonic - Boolean that indicates whether the flow is transonic (default FALSE)
        // nOuterCorr - Integer containing the number of outer correctors (default 1)
#       include "readPISOControls.H"
        // Evaluates mean and max Courant Number, and the velocity magnitude
        // Courant Number = U * DT / DX (Definition)
        // Evaluated as: (Flux / Area) * DT / Cell Distance
#       include "CourantNo.H"
            
#       include "readMeasData.H"

        // Pressure-velocity PISO corrector
        {
            // Momentum predictor
            // - LHS of the velocity equation (Left Hand Side):
            // - 1) Time derivative of velocity (implicit evaluation)
            // - 2) Convective term (divergence of the velocity flux times the velocity)
            // - 3) Divergence of the deviatoric part of the stress tensor (divergence of the shear rate tensor)
            // - Laplacian term is included in the divergence term
            fvVectorMatrix UEqn
            (
                fvm::ddt(U)
              + fvm::div(phi, U)
              + turbulence->divDevReff(U)
            );
            
            // - Relax the UEqn matrix according to the relaxation factor for U specified in controlDict
            UEqn.relax();
            
            // - Check if the momentum predictor (in fvSolution) is set to true
            if (momentumPredictor)
            {
                // - LHS of the velocity equation is set equal to the RHS side (Right Hand Side)
                // - Solve function uses the specified solver in fvSolution to find U (using the last known value of P)
                // - The found velocity field is NOT divergence free
                solve(UEqn == -fvc::grad(p));
            }
            
#           include "KalmanPredictionStep.H"
        
            // --- PISO loop (Inner Loop over specified number of correctors)

            for (int corr = 0; corr < nCorr; corr++)
            {
               
                // - UEqn.A returns the diagonal elements of the coefficient matrix
                // - A is the diagonal matrix of UEqn divided by the mesh single cell volume
                // - Note: A is function of U due to non-linearity of convection term               
                volScalarField rUA = 1.0/UEqn.A();

#               include "KalmanCorrectionStep.H"
                
                
                // - U is interpolated from cell centres to faces and then dotted with the face normals
                // - The second term accounts for the divergence of the face velocity field (stabilizing term)
                phi = (fvc::interpolate(U) & mesh.Sf())
                    + fvc::ddtPhiCorr(rUA, U, phi);

                // - Force inlet and outlet fluxes to obey continuity (necessary to find a solution for pressure)
                adjustPhi(phi, U, p);

               for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
                {
                    // Pressure corrector

                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian(rUA, p) == fvc::div(phi)
                    );

                    pEqn.setReference(pRefCell, pRefValue);

                    if
                    (
                        corr == nCorr-1
                     && nonOrth == nNonOrthCorr
                    )
                    {
                        pEqn.solve(mesh.solutionDict().solver("pFinal"));
                    }
                    else
                    {
                        pEqn.solve();
                    }

                    if (nonOrth == nNonOrthCorr)
                    {
                        phi -= pEqn.flux();
                    }
                }
                
                
                // - Calculates and print the continuity errors declared in InitContinuityErrs.H
#               include "continuityErrs.H"

#               include "divFreeCondition.H"
                
                U.correctBoundaryConditions();
                
#               include "calculateMisfit.H"
            }
        }
                
#       include "KalmanRegularisationStep.H"
        
        turbulence->correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
