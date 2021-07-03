
function saveResults2VTKv(
		filename::String, 
		testMesh::mesh2d, 
		testFields::fields2d)

	vtkfile = vtk_grid(filename, testMesh.xNodes, testMesh.yNodes, testMesh.VTKCells);
	vtk_point_data(vtkfile, testFields.densityNodes, "density");
	vtk_point_data(vtkfile, testFields.UxNodes, "Ux");
	vtk_point_data(vtkfile, testFields.UyNodes, "Uy");
	vtk_point_data(vtkfile, testFields.pressureNodes, "pressure");
	outfiles = vtk_save(vtkfile);
	
end


function createFields2d(testMesh::mesh2d, thermo::THERMOPHYSICS)


	##display("set initial and boundary conditions ...");

	phsLeftBC = zeros(Float64,4);
	phsTopBC = zeros(Float64,4);

	phsLeftBC[1] = 1.0;
	phsLeftBC[2] = 290.0;
	phsLeftBC[3] = 0.0;
	phsLeftBC[4] = 7143.0;

	phsTopBC[1] = 1.7;
	phsTopBC[2] = 263.72;
	phsTopBC[3] = -51.62;
	phsTopBC[4] = 15282.0;

	# solution intialization 
	# UphysCells = zeros(Float64, testMesh.nCells,4);
	#UconsCellsOld = zeros(Float64,testMesh.nCells,4);
	#UconsCellsNew = zeros(Float64,testMesh.nCells,4);

	densityCells = zeros(Float64, testMesh.nCells); 
	UxCells = zeros(Float64, testMesh.nCells); 
	UyCells = zeros(Float64, testMesh.nCells); 
	pressureCells = zeros(Float64, testMesh.nCells); 
	aSoundCells = zeros(Float64, testMesh.nCells); #speed of sound
	VMAXCells = zeros(Float64, testMesh.nCells); #max speed in domain

	#entropyCell = zeros(Float64, testMesh.nCells); #entropy

	densityNodes = zeros(Float64, testMesh.nNodes); 
	UxNodes = zeros(Float64, testMesh.nNodes); 
	UyNodes = zeros(Float64, testMesh.nNodes); 
	pressureNodes = zeros(Float64, testMesh.nNodes); 

	for i=1:testMesh.nCells

		densityCells[i] = phsLeftBC[1];
		UxCells[i] = phsLeftBC[2];
		UyCells[i] = phsLeftBC[3];
		pressureCells[i] = phsLeftBC[4];
			
		aSoundCells[i] = sqrt( thermo.Gamma * pressureCells[i]/densityCells[i] );
		VMAXCells[i]  = sqrt( UxCells[i]*UxCells[i] + UyCells[i]*UyCells[i] ) + aSoundCells[i];
		#entropyCell[i] = UphysCells[i,1]/(thermo.Gamma-1.0)*log(UphysCells[i,4]/UphysCells[i,1]*thermo.Gamma);
				
	end

	
	
	densityNodes = cells2nodesSolutionReconstructionWithStencils(testMesh, densityCells); 
	UxNodes = cells2nodesSolutionReconstructionWithStencils(testMesh, UxCells); 
	UyNodes = cells2nodesSolutionReconstructionWithStencils(testMesh, UyCells); 
	pressureNodes = cells2nodesSolutionReconstructionWithStencils(testMesh, pressureCells); 
		
	if (output.saveDataToVTK == 1)	
		filename = string("zzz",dynControls.curIter+1000);
		saveResults2VTK(filename, testMesh, densityNodes, "density");
		
	end

		


	# create fields 
	testFields2d = fields2d(
		densityCells,
		UxCells,
		UyCells,
		pressureCells,
		aSoundCells,
		VMAXCells,
		densityNodes,
		UxNodes,
		UyNodes,
		pressureNodes
		#UconsCellsOld,
		#UconsCellsNew
	);

	return testFields2d; 


end

function updateResidual!(
	Delta::Array{Float64,2},
	residualsVector1::Array{Any,1},
	residualsVector2::Array{Any,1},
	residualsVector3::Array{Any,1},
	residualsVector4::Array{Any,1},
	residualsVectorMax::Array{Float64,1},
	convergenceCriteria::Array{Float64,1},
	dynControls::DYNAMICCONTROLS,
	)

	residuals1::Float64   =sum( Delta[:,1].*Delta[:,1] );
	residuals2::Float64  = sum( Delta[:,2].*Delta[:,2] );
	residuals3::Float64  = sum( Delta[:,3].*Delta[:,3] );
	residuals4::Float64  = sum( Delta[:,4].*Delta[:,4] );
	
	push!(residualsVector1, residuals1);
	push!(residualsVector2, residuals2);
	push!(residualsVector3, residuals3);
	push!(residualsVector4, residuals4);
	
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



end


function updateVariablesM2(
	Delta::Array{Float64,2},
	UconsCellsOld::Array{Float64,2},
	UconsCellsNew::Array{Float64,2},
	testMesh::mesh2d,
	testfields2d::fields2d,
	dynControls::DYNAMICCONTROLS)
	
	for i=1:testMesh.nCells
	
	
		Delta[i,1] = UconsCellsNew[i,1] - UconsCellsOld[i,1];
		Delta[i,2] = UconsCellsNew[i,2] - UconsCellsOld[i,2];
		Delta[i,3] = UconsCellsNew[i,3] - UconsCellsOld[i,3];
		Delta[i,4] = UconsCellsNew[i,4] - UconsCellsOld[i,4];
	
		testfields2d.densityCells[i] = UconsCellsNew[i,1];
		testfields2d.UxCells[i] 	 = UconsCellsNew[i,2]/UconsCellsNew[i,1];
		testfields2d.UyCells[i] 	 = UconsCellsNew[i,3]/UconsCellsNew[i,1];
		testfields2d.pressureCells[i] = (thermo.Gamma-1.0)*( UconsCellsNew[i,4] - 0.5*( UconsCellsNew[i,2]*UconsCellsNew[i,2] + UconsCellsNew[i,3]*UconsCellsNew[i,3] )/UconsCellsNew[i,1] );

		testfields2d.aSoundCells[i] = sqrt( thermo.Gamma * testfields2d.pressureCells[i]/testfields2d.densityCells[i] );
		testfields2d.VMAXCells[i]  = sqrt( testfields2d.UxCells[i]*testfields2d.UxCells[i] + testfields2d.UyCells[i]*testfields2d.UyCells[i] ) + testfields2d.aSoundCells[i];
			
		UconsCellsOld[i,1] = UconsCellsNew[i,1];
		UconsCellsOld[i,2] = UconsCellsNew[i,2];
		UconsCellsOld[i,3] = UconsCellsNew[i,3];
		UconsCellsOld[i,4] = UconsCellsNew[i,4];
		
	end
	
	cells2nodesSolutionReconstructionWithStencilsImplicit!(testMesh, testfields2d); 
	
	(dynControls.rhoMax,id) = findmax(testfields2d.densityCells);
	(dynControls.rhoMin,id) = findmin(testfields2d.densityCells);
	
end

# DEPRICATED!!!
# function updateVariables!(
	# UconsCellsNew::Array{Float64,2},
	# testMesh::mesh2d,
	# testfields2d::fields2d,
	# dynControls::DYNAMICCONTROLS)
	
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
	
# end

function updateOutput!(
	timeVector::Array{Any,1},
	residualsVector1::Array{Any,1},
	residualsVector2::Array{Any,1},
	residualsVector3::Array{Any,1},
	residualsVector4::Array{Any,1},
	residualsVectorMax::Array{Float64,1}, 
	testMesh::mesh2d,
	testFields::fields2d,
	solControls::CONTROLS,
	output::outputCONTROLS,
	dynControls::DYNAMICCONTROLS)

	if (dynControls.verIter == output.verbosity)


		densityWarn = @sprintf("Density Min/Max: %f/%f", dynControls.rhoMin, dynControls.rhoMax);
		out = @sprintf("%0.6f\t %0.6f \t %0.6f \t %0.6f \t %0.6f \t %0.6f \t %0.6f", 
			dynControls.flowTime,
			dynControls.tau,
			residualsVector1[dynControls.curIter]./residualsVectorMax[1],
			residualsVector2[dynControls.curIter]./residualsVectorMax[2],
			residualsVector3[dynControls.curIter]./residualsVectorMax[3],
			residualsVector4[dynControls.curIter]./residualsVectorMax[4],
			dynControls.cpuTime
			 );
		#outputS = string(output, cpuTime);
		#println(outputS); 
		println(out); 
		println(densityWarn);
		
		#densityF = cells2nodesSolutionReconstructionWithStencils(testMesh,  testFields.densityCells  )
		
		if (output.saveDataToVTK == 1)
			filename = string("zzz",dynControls.curIter+1000); 
			
			
			saveResults2VTK(filename, testMesh, testFields.densityNodes, "density");
			
		end
		 
		thetaRad::Float64 = 61.0*pi/180.0;
		SL::Float64 = tan(thetaRad)*1.0;
		
		exactSolutionX = [
		0.0
		SL
		4.1];
		exactSolutionY = [
		1.0
		0.0
		1.0
		];
		
	
		if (solControls.plotResidual == 1)	


			
			subplot(2,1,1);
			
			cla();
			
			
			tricontourf(testMesh.xNodes,testMesh.yNodes, testFields.densityNodes,pControls.nContours,vmin=pControls.rhoMINcont,vmax=pControls.rhoMAXcont);
			#tricontour(testMesh.xNodes,testMesh.yNodes, testFields.densityNodes,pControls.nContours,vmin=pControls.rhoMINcont,vmax=pControls.rhoMAXcont);
			
			plot(exactSolutionX,exactSolutionY,"--k",linewidth = 1.5);
			
			set_cmap("jet");
			xlabel("x");
			ylabel("y");
			title("Contours of density");
			axis("equal");
			
			subplot(2,1,2);
			cla();
			#GR.clearws();
			#plot(timeVector, residualsVector1./residualsVectorMax[1], timeVector, residualsVector2./residualsVectorMax[2],timeVector,residualsVector3./	residualsVectorMax[3],timeVector, residualsVector4./residualsVectorMax[4]); 
			plot(timeVector, residualsVector1./residualsVectorMax[1],"-r",label="continuity"); 
			plot(timeVector, residualsVector2./residualsVectorMax[2],"-g",label="momentum ux"); 
			plot(timeVector, residualsVector3./residualsVectorMax[3],"-b",label="momentum uy"); 
			plot(timeVector, residualsVector4./residualsVectorMax[4],"-c",label="energy"); 
			
			yscale("log");	
			xlabel("flow time [s]");
			ylabel("Res");
			title("Residuals");
			legend();
			
			pause(1.0e-5);
			
			
		end

   

		#pause(1.0e-3);
		dynControls.verIter = 0; 

	end


end

