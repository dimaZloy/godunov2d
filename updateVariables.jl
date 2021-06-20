#display("done ...");


if (nprocs() == 1) #serial 

	
	Delta =  deepcopy(UconsCellsNew - UconsCellsOld); 
	UphysCells = deepcopy(cns2dphs2d(UconsCellsNew,testMesh.nCells,thermo.Gamma));
	#UphysNodes = deepcopy( cells2nodesSolutionReconstructionWithStencilsM(nCells, nNodes, nNeibCells, UphysCells , cell_clusters, node_stencils)); 
	UphysNodes = deepcopy( cells2nodesSolutionReconstructionWithStencils(testMesh, UphysCells )  ); 
 
	#Delta =  UconsCellsNew - UconsCellsOld; 
	#UphysCells = cns2dphs2d(UconsCellsNew,thermo.Gamma);
	#UphysNodes = cells2nodesSolutionReconstructionWithStencilsM(nCells, nNodes, nNeibCells, UphysCells , cell_clusters, node_stencils); 

	for i=1:testMesh.nCells
    	aSoundCell[i] = sqrt( thermo.Gamma * UphysCells[i,4]/UphysCells[i,1] );
      	VMAXCell[i]  = sqrt( UphysCells[i,2]*UphysCells[i,2] + UphysCells[i,3]*UphysCells[i,3] ) + aSoundCell[i];
      	densityCell[i] = UphysCells[i,1];
		entropyCell[i] = UphysCells[i,1]/(thermo.Gamma-1.0)*log(UphysCells[i,4]/UphysCells[i,1]*thermo.Gamma);
	end


	#densityNodes = cells2nodesSolutionReconstructionWithStencilsV(nCells, nNodes, nNeibCells, density , cell_clusters, node_stencils); 
	densityNodes = cells2nodesSolutionReconstructionWithStencils(testMesh, densityCell); 

	(dynControls.rhoMax,id) = findmax(densityCell);
	(dynControls.rhoMin,id) = findmin(densityCell);

	UconsCellsOld =  deepcopy(UconsCellsNew);
else

	

	# UphysCellsX = cns2dphs2dX(UconsCellsNewX,nCells,thermo.Gamma);
	# UphysNodesX = cells2nodesSolutionReconstructionWithStencilsMX(nCells, nNodes, nNeibCells, UphysCellsX , cell_clusters, node_stencils); 
 

	# for i=1:nCells
    	# aSoundX[i] = sqrt( thermoX.Gamma * UphysCellsX[i,4]/UphysCellsX[i,1] );
      	# VMAXX[i]  = sqrt( UphysCellsX[i,2]*UphysCellsX[i,2] + UphysCellsX[i,3]*UphysCellsX[i,3] ) + aSoundX[i];
      	# densityX[i] = UphysCellsX[i,1];
		# entropyX[i] = UphysCellsX[i,1]/(thermo.Gamma-1.0)*log(UphysCellsX[i,4]/UphysCellsX[i,1]*thermoX.Gamma);
	# end


	# densityNodesX = cells2nodesSolutionReconstructionWithStencilsV(nCells, nNodes, nNeibCells, densityX , cell_clusters, node_stencils); 
	# entropyNodesX = cells2nodesSolutionReconstructionWithStencilsV(nCells, nNodes, nNeibCells, entropyX , cell_clusters, node_stencils); 

	# (dynControls.rhoMax,id) = findmax(densityX);
	# (dynControls.rhoMin,id) = findmin(densityX);


end

