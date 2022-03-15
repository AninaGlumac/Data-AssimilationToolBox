                 // - Kalman correction (measurement available)               
                if (exists(runTime.timeName()/"UMeasure"))
                {
//U_test1=U;
                    U = UEqn().H()/AU;
//U_test2=U;
                    U.internalField() = U.internalField() +KalmanGain*(UField[0].internalField()-(HTrans*U.internalField()));
		
		    //U = U +KalmanGain*(UField[0]-(HTrans*U));
//U_test3=U;
                }
                // - PISO loop (no measurement available)
                else
                {
                    U = UEqn().H()/AU;
                }
