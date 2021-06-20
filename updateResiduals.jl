

if (nprocs() == 1) #serial 
	residuals1   =sum( Delta[:,1].*Delta[:,1] );
	residuals2  = sum( Delta[:,2].*Delta[:,2] );
	residuals3  = sum( Delta[:,3].*Delta[:,3] );
	residuals4  = sum( Delta[:,4].*Delta[:,4] );
else # parallel
	residuals1   =sum( DeltaX[:,1].*DeltaX[:,1] );
	residuals2  = sum( DeltaX[:,2].*DeltaX[:,2] );
	residuals3  = sum( DeltaX[:,3].*DeltaX[:,3] );
	residuals4  = sum( DeltaX[:,4].*DeltaX[:,4] );	 
end

	residualsVector1 = [residualsVector1; residuals1];
	residualsVector2 = [residualsVector2; residuals2];
	residualsVector3 = [residualsVector3; residuals3];
	residualsVector4 = [residualsVector4; residuals4];

	if (dynControls.curIter<6 && dynControls.curIter>1)

   		(residualsVectorMax[1],id1) = findmax(residualsVector1[1:dynControls.curIter]);	
   		(residualsVectorMax[2],id2) = findmax(residualsVector2[1:dynControls.curIter]);		
   		(residualsVectorMax[3],id3) = findmax(residualsVector3[1:dynControls.curIter]);		
   		(residualsVectorMax[4],id4) = findmax(residualsVector4[1:dynControls.curIter]);		

	end

	if ( (dynControls.curIter>5) && 
    	(residualsVector1[dynControls.curIter]./residualsVectorMax[1] <= convergenceCriteria[1]) &&
     	(residualsVector2[dynControls.curIter]./residualsVectorMax[2] <= convergenceCriteria[2]) &&
     	(residualsVector3[dynControls.curIter]./residualsVectorMax[3] <= convergenceCriteria[3]) &&
     	(residualsVector4[dynControls.curIter]./residualsVectorMax[4] <= convergenceCriteria[4]) )

	 	dynControls.isSolutionConverged  = 1; 

	end

	


if (dynControls.flowTime >= solControls.stopTime)
	dynControls.isRunSimulation = 0;
end