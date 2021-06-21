

# ## DEPRICATED!!!

# timeVector = [];
# residualsVector1 = []; 
# residualsVector2 = []; 
# residualsVector3 = []; 
# residualsVector4 = []; 

# Delta2 = zeros(Float64,testMesh.nCells,4); 
# residualsVectorMax = ones(Float64,4);
# convergenceCriteria = [1e-5;1e-5;1e-5;1e-5;];


# display("set initial and boundary conditions ...");

# phsLeftBC = zeros(Float64,4);
# phsTopBC = zeros(Float64,4);

# phsLeftBC[1] = 1.0;
# phsLeftBC[2] = 290.0;
# phsLeftBC[3] = 0.0;
# phsLeftBC[4] = 7143.0;

# phsTopBC[1] = 1.7;
# phsTopBC[2] = 263.72;
# phsTopBC[3] = -51.62;
# phsTopBC[4] = 15282.0;

# # solution intialization 
# UphysCells = zeros(Float64, testMesh.nCells,4);
# UconsCellsOld = zeros(Float64,testMesh.nCells,4);
# aSoundCell = zeros(Float64, testMesh.nCells); #speed of sound
# VMAXCell = zeros(Float64, testMesh.nCells); #max speed in domain
# densityCell = zeros(Float64, testMesh.nCells); #density 
# entropyCell = zeros(Float64, testMesh.nCells); #entropy


# for i=1:testMesh.nCells
	# UphysCells[i,1] = phsLeftBC[1];
	# UphysCells[i,2] = phsLeftBC[2];
	# UphysCells[i,3] = phsLeftBC[3];
	# UphysCells[i,4] = phsLeftBC[4];
		
    # aSoundCell[i] = sqrt( thermo.Gamma * UphysCells[i,4]/UphysCells[i,1] );
    # VMAXCell[i]  = sqrt( UphysCells[i,2]*UphysCells[i,2] + UphysCells[i,3]*UphysCells[i,3] ) + aSoundCell[i];
    # densityCell[i] = UphysCells[i,1];
	# entropyCell[i] = UphysCells[i,1]/(thermo.Gamma-1.0)*log(UphysCells[i,4]/UphysCells[i,1]*thermo.Gamma);
			
# end

# UphysNodes = cells2nodesSolutionReconstructionWithStencils(testMesh, UphysCells ); 

# densityNodes = zeros(Float64, testMesh.nNodes);
# densityNodes = UphysNodes[:,1]; 
	
# if (output.saveDataToVTK == 1)	
	# filename = string("zzz",dynControls.curIter+1000);
	# saveResults2VTK(filename, testMesh, densityNodes, "density");
# end

	
# UconsCellsOld = phs2dcns2d(UphysCells,thermo.Gamma); #old vector

	
# UconsDelta = deepcopy(UconsCellsOld); #delta - residual vector
# UconsCellsN1 = deepcopy(UconsCellsOld); #RK4 i-stage
# UconsCellsN2 = deepcopy(UconsCellsOld); #RK4 i-stage
# UconsCellsN3 = deepcopy(UconsCellsOld); #RK4 i-stage
# UconsCellsNew = deepcopy(UconsCellsOld); #new  vector 


