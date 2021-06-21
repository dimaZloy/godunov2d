
# DEPRICATED!!!!

# if (solControls.timeStepMethod == 1)
	# dynControls.flowTime += dynControls.tau;  
# else
	# dynControls.flowTime += solControls.dt;  
# end
# timeVector = [timeVector; dynControls.flowTime]; 

# dynControls.curIter += 1; 
# dynControls.verIter += 1;


# include("updateVariables.jl")
# include("updateResiduals.jl");
# include("updateOutput.jl");


# if (flowTime>= solControls.stopTime || dynControls.isSolutionConverged == 1)
	# dynControls.isRunSimulation = 0;
	
	# if (dynControls.isSolutionConverged == true)
		# println("Solution converged! ");
	# else
		# println("Simultaion flow time reached the set Time!");
	# end
	
	# if (output.saveResiduals == 1)
		# println("Saving Residuals ... ");
		# cd(dynControls.localTestPath);
		# saveResiduals(output.fileNameResiduals, timeVector, residualsVector1, residualsVector2, residualsVector3, residualsVector4);
		# cd(dynControls.globalPath);
	# end
	# if (output.saveResults == 1)
		# println("Saving Results ... ");
		# cd(dynControls.localTestPath);
		# saveSolution(output.fileNameResults, testMesh.xNodes, testMesh.yNodes, UphysNodes);
		# cd(dynControls.globalPath);
	# end
	
# end


