
#function SecondOrderUpwindM2(bettaKJ::Float64, UConsCellsOld::Array{Float64,2}, 
#	testMesh::mesh2d, testFields::fields2d, thermo::THERMOPHYSICS, solControls::CONTROLS, dynControls::DYNAMICCONTROLS, UconsCellsNew::Array{Float64,2})
	
function SecondOrderUpwindM2(bettaKJ::Float64, dt::Float64, UConsCellsOld::Array{Float64,2}, 
	testMesh::mesh2d, testFields::fields2d, thermo::THERMOPHYSICS, UconsCellsNew::Array{Float64,2})

	
	uLeftp = zeros(Float64,4);
	FLUXES = zeros(Float64,testMesh.nCells,4);
	
	for i=1:testMesh.nCells
    

	
		num_nodes::Int64 = testMesh.mesh_connectivity[i,3];
		
		
		uLeftp[1] = testFields.densityCells[i];
		uLeftp[2] = testFields.UxCells[i];
		uLeftp[3] = testFields.UyCells[i];
		uLeftp[4] = testFields.pressureCells[i];
		
		# uRightp[1] = testFields.densityCells[i];
		# uRightp[2] = testFields.UxCells[i];
		# uRightp[3] = testFields.UyCells[i];
		# uRightp[4] = testFields.pressureCells[i];
	   
		
		if (num_nodes == 3)
		
			edge_flux1 = ( computeInterfaceSlope(i, 1, testMesh, testFields, uLeftp) );
			edge_flux2 = ( computeInterfaceSlope(i, 2, testMesh, testFields, uLeftp) );
			edge_flux3 = ( computeInterfaceSlope(i, 3, testMesh, testFields, uLeftp) );
				
			FLUXES[i,1] = edge_flux1[1] + edge_flux2[1] + edge_flux3[1];
			FLUXES[i,2] = edge_flux1[2] + edge_flux2[2] + edge_flux3[2];
			FLUXES[i,3] = edge_flux1[3] + edge_flux2[3] + edge_flux3[3];
			FLUXES[i,4] = edge_flux1[4] + edge_flux2[4] + edge_flux3[4];
			

		elseif (num_nodes == 4)
			
			edge_flux1 = ( computeInterfaceSlope(i, 1, testMesh, testFields, uLeftp) );
			edge_flux2 = ( computeInterfaceSlope(i, 2, testMesh, testFields, uLeftp) );
			edge_flux3 = ( computeInterfaceSlope(i, 3, testMesh, testFields, uLeftp) );
			edge_flux4 = ( computeInterfaceSlope(i, 4, testMesh, testFields, uLeftp) );
				
			FLUXES[i,1] = edge_flux1[1] + edge_flux2[1] + edge_flux3[1] + edge_flux4[1];
			FLUXES[i,2] = edge_flux1[2] + edge_flux2[2] + edge_flux3[2] + edge_flux4[2];
			FLUXES[i,3] = edge_flux1[3] + edge_flux2[3] + edge_flux3[3] + edge_flux4[3];
			FLUXES[i,4] = edge_flux1[4] + edge_flux2[4] + edge_flux3[4] + edge_flux4[4];

		else
			
			display("something wrong in flux calculations ... ")
			
		end
		
		
		UconsCellsNew[i,1] = ( UConsCellsOld[i,1] - FLUXES[i,1]*bettaKJ*dt*testMesh.Z[i] );
		UconsCellsNew[i,2] = ( UConsCellsOld[i,2] - FLUXES[i,2]*bettaKJ*dt*testMesh.Z[i] );
		UconsCellsNew[i,3] = ( UConsCellsOld[i,3] - FLUXES[i,3]*bettaKJ*dt*testMesh.Z[i] );
		UconsCellsNew[i,4] = ( UConsCellsOld[i,4] - FLUXES[i,4]*bettaKJ*dt*testMesh.Z[i] );
   
		# if (solControls.timeStepMethod == 1)
		
			# UconsCellsNew[i,1] = ( UConsCellsOld[i,1] - FLUXES[i,1]*bettaKJ*testMesh.Z[i]*dynControls.tau );
			# UconsCellsNew[i,2] = ( UConsCellsOld[i,2] - FLUXES[i,2]*bettaKJ*testMesh.Z[i]*dynControls.tau );
			# UconsCellsNew[i,3] = ( UConsCellsOld[i,3] - FLUXES[i,3]*bettaKJ*testMesh.Z[i]*dynControls.tau );
			# UconsCellsNew[i,4] = ( UConsCellsOld[i,4] - FLUXES[i,4]*bettaKJ*testMesh.Z[i]*dynControls.tau );
			
		# else
		
			# UconsCellsNew[i,1] = ( UConsCellsOld[i,1] - FLUXES[i,1]*bettaKJ*testMesh.Z[i]*solControls.dt );
			# UconsCellsNew[i,2] = ( UConsCellsOld[i,2] - FLUXES[i,2]*bettaKJ*testMesh.Z[i]*solControls.dt );
			# UconsCellsNew[i,3] = ( UConsCellsOld[i,3] - FLUXES[i,3]*bettaKJ*testMesh.Z[i]*solControls.dt );
			# UconsCellsNew[i,4] = ( UConsCellsOld[i,4] - FLUXES[i,4]*bettaKJ*testMesh.Z[i]*solControls.dt );
			
		# end
   
   
	end # i - loop for all cells

	
	##return UconsCellsNew;

end


