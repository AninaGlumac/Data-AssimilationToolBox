            volScalarField HMatrix = 1.0/(UEqn.A()*runTime.deltaTValue());

            volScalarField TurbulentC
            (
                "TurbulentC",
                turbulence->nut()/turbulence->nu()
            );
            

            //  - Model constant for the turbulence case (depends on C)
            // - QModel = Model variance matrix (supposed diagonal)
            // - PCovariance = Covariance matrix

            TurbulentC = 1+TurbulentC;
            QModel = C*TurbulentC*runTime.deltaTValue()*deltaX;
            PCovariance = HMatrix*PCovarianceOld*HMatrix*pomDim+QModel;

            SMeasure = HTrans*PCovariance*HTrans+MeasureVariance;
        
            // - Kalman Gain evaluated only if measurement exists
            if (exists(runTime.timeName()/"UMeasure"))
            {
                    // - Kalman Gain = measurement importance weight
                    KalmanGain = PCovariance*HTrans/SMeasure;
                    
                    //Tried to multiplicate with Flag field in order to set 0 to all possitions that not related to observation points, but not working, so it is done through the loop
                    forAll(mesh.C(), i)
                    {
                        if (Flag.internalField()[i]==0)
                        {
                            SMeasure.internalField()[i]=0;
                            HTrans.internalField()[i]=0;
                            KalmanGain.internalField()[i]=0;
                        }
                     }//KalmanGain=KalmanGain*Flag;
                     
                    KSum=gSum(KalmanGain);//*Flag);
                    autoPtr<OFstream> file_Kalman;

//                     if (Pstream::master())
//                     {
//                         file_Kalman.reset( new OFstream(runTime.path()/"KalmanGainN"));
//                     }
                    
                    if (Pstream::master())
                    {
                        Info<<runTime.time().value()<<" "<<KSum<<endl;
                        
//                         file_Kalman()<<KSum/NMeasure<< endl;
                        OFstream strKA(runTime.path()/"KalmanGain", ios_base::app);
                        
                        strKA << (KSum/NMeasure) << endl;
                    }
            }
