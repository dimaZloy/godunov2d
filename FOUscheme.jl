


	
function FirstOrderUpwindM2(bettaKJ::Float64, dt::Float64, UConsCellsOld::Array{Float64,2}, 
	testMesh::mesh2d, testFields::fields2d, thermo::THERMOPHYSICS, UconsCellsNew::Array{Float64,2})
	

	
	uLeftp = zeros(Float64,4);
	uRightp = zeros(Float64,4);	

	uConsLeftp = zeros(Float64,4);
	uConsRightp = zeros(Float64,4);	

	FLUXES = zeros(Float64,testMesh.nCells,4);
	
	for i = 1:testMesh.nCells
	
		ck::Int64 = testMesh.mesh_connectivity[i,3]; 

		uLeftp[1] = testFields.densityCells[i];
		uLeftp[2] = testFields.UxCells[i];
		uLeftp[3] = testFields.UyCells[i];
		uLeftp[4] = testFields.pressureCells[i];
		
		uRightp[1] = testFields.densityCells[i];
		uRightp[2] = testFields.UxCells[i];	
		uRightp[3] = testFields.UyCells[i];	
		uRightp[4] = testFields.pressureCells[i];
		
		
		
		for k =1:ck
					
			side::Float64 = testMesh.cell_edges_length[i,k];
			nx::Float64   = testMesh.cell_edges_Nx[i,k];
			ny::Float64   = testMesh.cell_edges_Ny[i,k];

			ek::Int64 = testMesh.cell_stiffness[i,k]; 
			edge_flux = zeros(Float64,4);


			if (ek>=1 && ek<=testMesh.nCells ) 
				
				uRightp[1] = testFields.densityCells[ek];
				uRightp[2] = testFields.UxCells[ek];
				uRightp[3] = testFields.UyCells[ek];
				uRightp[4] = testFields.pressureCells[ek];
				
				
			else
				uRightp = ComputeUPhysFromBoundaries(i,k, ek, uRightp, nx,ny);
				
			end 
			
			# flux approximation:
			#  0 - Roe
			#  1 - exact Riemann solver
			#  2 - AUSM+ 
			#  3 - AUSM+up
			
			
			#edge_flux = RoeFlux2d(uRightp,uLeftp, nx,ny,side,thermo.Gamma);
			edge_flux = AUSMplusFlux2d(uRightp,uLeftp,nx,ny,side,thermo.Gamma);
			

			#FLUXES[i,:] += edge_flux;
			FLUXES[i,1] = FLUXES[i,1] + edge_flux[1];
			FLUXES[i,2] = FLUXES[i,2] + edge_flux[2];
			FLUXES[i,3] = FLUXES[i,3] + edge_flux[3];
			FLUXES[i,4] = FLUXES[i,4] + edge_flux[4];
			
			
		end # K for neib cells 
			
		UconsCellsNew[i,1] = ( UConsCellsOld[i,1] - FLUXES[i,1]*bettaKJ*dt*testMesh.Z[i] );
		UconsCellsNew[i,2] = ( UConsCellsOld[i,2] - FLUXES[i,2]*bettaKJ*dt*testMesh.Z[i] );
		UconsCellsNew[i,3] = ( UConsCellsOld[i,3] - FLUXES[i,3]*bettaKJ*dt*testMesh.Z[i] );
		UconsCellsNew[i,4] = ( UConsCellsOld[i,4] - FLUXES[i,4]*bettaKJ*dt*testMesh.Z[i] );
		
	end	#  i for loop cells 

	
end