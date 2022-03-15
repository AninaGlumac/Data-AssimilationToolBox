        if (exists(runTime.timeName()/"UMeasure"))
        {
            // - Kalman regularization step (only if measurement available)
            
            // - The covariance matrix is updated with the Kalman Gain
            PCovariance = (1 - KalmanGain*HTrans)*PCovarianceOld;
            PCovarianceOld = PCovariance;
        }
         // - L2 Error Evaluation
        scalar sumVelocityError = fvc::domainIntegrate((Misfit)&(Misfit)).value();
        scalar VelocityMeasure = 0;

        if (exists(runTime.timeName()/"UMeasure"))
        {
            VelocityMeasure = fvc::domainIntegrate(UField[0]&UField[0]).value();
            Info<< "Absolute Velocity Field Error = " << Foam::sqrt(sumVelocityError) << endl;
            Info << "Relative Velocity Field Error = " << Foam::sqrt(sumVelocityError/VelocityMeasure) << endl;
            Info << " Velocity Field Measure = " << Foam::sqrt(VelocityMeasure) << endl;
            OFstream strL(runTime.path()/"VelocityL2Error", ios_base::app);
            strL << Foam::sqrt(sumVelocityError/VelocityMeasure) << endl;
        }
