        // - UField is initialized with the current velocity, then it is set equal to zero
        UField.set(0, new volVectorField(U, mesh));
        UField[0].internalField() = UTmp[0].internalField();
        
          // - Reads the measure (only if present), which in this case is written in the standard OF format
        if (exists(runTime.timeName()/"UMeasure"))
        {
            NMeasure = 0;
            IOobject UReader
            (
                "UMeasure",
                runTime.timeName(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            );
        
            UField.set(0, new volVectorField(UReader,mesh));
            forAll(mesh.C(), i)
            {
                if (UField[0].internalField()[i].component(1) != 0)
                {
                    // - Flag is equal to 1 in cells where the measurement is present
                    Flag.internalField()[i] = 1;
                    NMeasure=NMeasure + 1;
                }
            }
            
            reduce(NMeasure, sumOp<scalar>());
            
            if (Pstream::master())
            {
                Info<<runTime.time().value()<<" "<<NMeasure<<endl;
            }
            Misfit = UField[0];
        }
