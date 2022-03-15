    //open file couplingProperties in directory: casename/constant
    IOdictionary DAProperties
    (
        IOobject
        (
            "DAProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );



    //get model constants
    scalar C(readScalar(DAProperties.lookup("C_Model")));
    //scalar os(readScalar(DAProperties.lookup("os_Model"))); - not taken into account, taken =1 so it is on a conservative side
    //scalar ot(readScalar(DAProperties.lookup("os_Model")));

