                 // - Kalman correction (measurement available)               
                if (exists(runTime.timeName()/"UMeasure"))
                {
                    
                    U.internalField() = U.internalField() +KalmanGain*(UField[0].internalField()-(HTrans*U.internalField()));
                    
//UBACILA OVO ISPOD!!!!
                    rUA = 1.0/UEqn.A();
                    
                }
                // - PISO loop (no measurement available)
                else
                {
                    U = rUA*UEqn.H();
                }
