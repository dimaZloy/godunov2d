
## DEPRICATED!!!

# #display("done ...");


# if (nprocs() == 1) #serial 

	
	# Delta =  deepcopy(UconsCellsNew - UconsCellsOld); 
	# #UphysCells = deepcopy(cns2dphs2d(UconsCellsNew,testMesh.nCells,thermo.Gamma));
	
	# #display(UconsCellsNew)
	
	# for i=1:testMesh.nCells


	
		# testfields2d.densityCells[i] = UconsCellsNew[i,1];
		# testfields2d.UxCells[i] 	 = UconsCellsNew[i,2]/UconsCellsNew[i,1];
		# testfields2d.UyCells[i] 	 = UconsCellsNew[i,3]/UconsCellsNew[i,1];
		# testfields2d.pressureCells[i] = (thermo.Gamma-1.0)*( UconsCellsNew[i,4] - 0.5*( UconsCellsNew[i,2]*UconsCellsNew[i,2] + UconsCellsNew[i,3]*UconsCellsNew[i,3] )/UconsCellsNew[i,1] );

		
		# testfields2d.aSoundCells[i] = sqrt( thermo.Gamma * testfields2d.pressureCells[i]/testfields2d.densityCells[i] );
		# testfields2d.VMAXCells[i]  = sqrt( testfields2d.UxCells[i]*testfields2d.UxCells[i] + testfields2d.UyCells[i]*testfields2d.UyCells[i] ) + testfields2d.aSoundCells[i];
		
	# end
	
	# cells2nodesSolutionReconstructionWithStencilsImplicit!(testMesh, testfields2d); 
	
	# (dynControls.rhoMax,id) = findmax(testfields2d.densityCells);
	# (dynControls.rhoMin,id) = findmin(testfields2d.densityCells);

	# UconsCellsOld =  deepcopy(UconsCellsNew);
# else

	

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

