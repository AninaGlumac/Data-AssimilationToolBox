    Info << "Reading field p\n" << endl;
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

    Info << "Reading field U\n" << endl;
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
volVectorField U_test1
    (
        IOobject
        (
            "U_test1",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("U_test1", dimensionSet(0,1,-1,0,0,0,0),vector::zero)
    );

volVectorField U_test2
    (
        IOobject
        (
            "U_test2",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("U_test2", dimensionSet(0,1,-1,0,0,0,0),vector::zero)
    );

volVectorField U_test3
    (
        IOobject
        (
            "U_test3",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("U_test3", dimensionSet(0,1,-1,0,0,0,0),vector::zero)
    );
#   include "createPhi.H"


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue);


    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::RASModel> turbulence
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );
    volScalarField PCovarianceOld
    (
        IOobject
        (
            "PCovarianceOld",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
            //IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("PCovarianceOld", dimensionSet(0,0,0,0,0,0,0),1)
    );
    
    volScalarField PCovariance
    (
        IOobject
        (
            "PCovariance",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
            //IOobject::AUTO_WRITE            
        ),
        mesh,
        dimensionedScalar("PCovariance", dimensionSet(0,0,0,0,0,0,0),1)
    );
    
    volScalarField QModel
    (
        IOobject
        (
            "QModel",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
            //IOobject::AUTO_WRITE           
        ),
        mesh,
        dimensionedScalar("QModel", dimensionSet(0,0,0,0,0,0,0),0)
    );
    
    volScalarField pomDim
    (
        IOobject
        (
            "pomDim",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
            //IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("pomDim", dimensionSet(0,0,-2,0,0,0,0),1)
    );
    
    volScalarField deltaX
    (
        IOobject
        (
            "deltaX",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
            //IOobject::AUTO_WRITE
        ),
        mesh,
      dimensionedScalar("deltaX", dimensionSet(0,0,0,0,0,0,0),1)
    );
    
    volScalarField MeasureVariance
    (
        IOobject
        (
            "MeasureVariance",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            //IOobject::NO_WRITE
            IOobject::AUTO_WRITE
        ),
        mesh       
    );
    
    volScalarField HTrans
    (
        IOobject
        (
            "HTrans",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
            //IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("HTrans", dimensionSet(0,0,0,0,0,0,0),1)
    );
    
    volScalarField SMeasure
    (
        IOobject
        (
            "SMeasure",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
            //IOobject::AUTO_WRITE
        ),
        mesh,
      dimensionedScalar("SMeasure", dimensionSet(0,0,0,0,0,0,0),0)
    );
      
    //volScalarField HMatrix
    //(
    //    mesh,
    //    dimensionedScalar("zero", p.dimensions(), 0)
    //)
    volScalarField KalmanGain
    (
        IOobject
        (
            "KalmanGain",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("KalmanGain", dimensionSet(0,0,0,0,0,0,0),0)
    );
    
    volScalarField Flag
    (
        IOobject
        (
            "Flag",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Flag", dimensionSet(0,0,0,0,0,0,0),0)
    );

    volVectorField Misfit
    (
        IOobject
        (            
            "Misfit",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("Misfit", dimensionSet(0,1,-1,0,0,0,0),vector::zero)
    );
