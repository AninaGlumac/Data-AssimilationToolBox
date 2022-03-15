                if (exists(runTime.timeName()/"UMeasure"))
                {
                    forAll(mesh.C(), i)
                    {
                        if (Flag.internalField()[i]==1)
                        {
                            Misfit.internalField()[i] = cmptMag(UField[0].internalField()[i] - U.internalField()[i]);
                        }
                    }
                    double maxMisX= max(Misfit.internalField().component(0));
                    
                    reduce(maxMisX, maxOp<scalar>());
                    
                    
                    double maxXmeas= max(UField[0].internalField().component(0));
                    
                    reduce(maxXmeas, maxOp<scalar>());
                    
                    double maxMisY= max(Misfit.internalField().component(1));
                    
                    reduce(maxMisY, maxOp<scalar>());
                    
                    
                    double maxYmeas= max(UField[0].internalField().component(1));
                    
                    reduce(maxYmeas, maxOp<scalar>());
                    
                    if (Pstream::master())
                    {
                        Info << "Maximum Misfit X = " << maxMisX << endl;
                        
                        OFstream strMX(runTime.path()/"MaxMisfitX", ios_base::app);
                        
                        strMX << maxMisX/maxXmeas << endl;
                        
                                                       
                        Info << "Maximum Misfit Y = " << maxMisX << endl;
                        
                        OFstream strMY(runTime.path()/"MaxMisfitY", ios_base::app);
			//Info<<"Ovde sam!"<<endl;
                        
                        //strMY << maxMisY/maxYmeas << endl;
    
                        //Info<<"Ovde sam!"<<endl;

                    }
                }
