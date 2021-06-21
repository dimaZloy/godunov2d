
# DEPRICATED!!!

# if (dynControls.verIter == output.verbosity)


	# densityWarn = @sprintf("Density Min/Max: %f/%f", dynControls.rhoMin, dynControls.rhoMax);
	# out = @sprintf("%0.6f\t %0.6f \t %0.6f \t %0.6f \t %0.6f \t %0.6f \t %0.6f", 
		# dynControls.flowTime,
		# dynControls.tau,
		# residualsVector1[dynControls.curIter]./residualsVectorMax[1],
		# residualsVector2[dynControls.curIter]./residualsVectorMax[2],
		# residualsVector3[dynControls.curIter]./residualsVectorMax[3],
		# residualsVector4[dynControls.curIter]./residualsVectorMax[4],
		# dynControls.cpuTime
		 # );
	# #outputS = string(output, cpuTime);
	# #println(outputS); 
	# println(out); 
	# println(densityWarn);
	
	
	# if (output.saveDataToVTK == 1)
		# filename = string("zzz",dynControls.curIter+1000); 
		# saveResults2VTK(filename, testMesh, testfields2d.densityNodes, "density");
	# end
	 
	
		

		
		
	
	# if (solControls.plotResidual == 1)	


		
		# subplot(2,1,1);
		
		# cla();
		
		# tricontourf(testMesh.xNodes,testMesh.yNodes,testfields2d.densityNodes,pControls.nContours,vmin=pControls.rhoMINcont,vmax=pControls.rhoMAXcont);
		
		# #if (pControls.contoursType == 1)
		# #	tricontourf(xNodes,yNodes,densityNodes,pControls.nContours);
		# #else
		# #	tricontour(xNodes,yNodes,densityNodes,pControls.nContours);
		# #end
		# #tricontour(xNodes,yNodes,densityNodes,pControls.nContours);
		# #clim(vmin=pControls.rhoMINcont,vmax=pControls.rhoMAXcont);	
		# #colorbar();
		# set_cmap("jet");
		# xlabel("x");
		# ylabel("y");
		# title("Contours of density");
		# axis("equal");
		
		# subplot(2,1,2);
		# cla();
		# #GR.clearws();
		# #plot(timeVector, residualsVector1./residualsVectorMax[1], timeVector, residualsVector2./residualsVectorMax[2],timeVector,residualsVector3./	residualsVectorMax[3],timeVector, residualsVector4./residualsVectorMax[4]); 
		# plot(timeVector, residualsVector1./residualsVectorMax[1],"-r",label="continuity"); 
		# plot(timeVector, residualsVector2./residualsVectorMax[2],"-g",label="momentum ux"); 
		# plot(timeVector, residualsVector3./residualsVectorMax[3],"-b",label="momentum uy"); 
		# plot(timeVector, residualsVector4./residualsVectorMax[4],"-c",label="energy"); 
		
		# yscale("log");	
		# xlabel("flow time [s]");
		# ylabel("Res");
		# title("Residuals");
		# legend();
		
		# pause(1.0e-5);
		
		
	# end

   

    # #pause(1.0e-3);
	# dynControls.verIter = 0; 

# end


