                //if (exists(runTime.timeName()/"UMeasure"))
                //{
                    // reinforces div free! (it really changes the U field!)
                //    U = UEqn.H()*rUA - rUA*fvc::grad(p);
                //}
                //else
                //{
                    U -= rUA*fvc::grad(p);
                //}
