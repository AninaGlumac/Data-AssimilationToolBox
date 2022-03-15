    
    
    // - UField is the container for the measurement
    // - UTmp is initialized with the initial velocity (which is zero)
    PtrList<volVectorField> UField(1);
    PtrList<volVectorField> UTmp(1);
    PtrList<volVectorField> USteady(1);

    UTmp.set(0,new volVectorField(U, mesh));

    // - NMeasure = number of locations the measure is present
    scalar NMeasure = 0;
    double KSum = 0;
    
    forAll (mesh.C(),i)
    {
        deltaX.internalField()[i]=mesh.V()[i]/mesh.magSf()[i];
    }
    
    // - C = User defined constant for the evaluation of Q (MUST BE OPTIMIZED!!!)
    const double C = 1;
