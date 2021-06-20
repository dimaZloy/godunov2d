


(dynControls.velmax,id) = findmax(VMAXCell);
dynControls.tau = solControls.CFL * testMesh.maxArea/(max(dynControls.velmax,1.0e-6));


UconsCellsNew = deepcopy(FirstOrderUpwindM1(1.0,UconsCellsOld, UphysCells ));

# if (solver.TimeDiscretization == 2) # 4th RK 

	# if (nprocs() == 1) #serial 
	
		# #UconsCellsN1 =  deepcopy( calculateVariablesAtStage(1.0/4.0,UconsCellsOld,UphysCells, UphysNodes) );
		# #UconsCellsN2 =  deepcopy( calculateVariablesAtStage(1.0/3.0,UconsCellsN1,UphysCells, UphysNodes) );
		# #UconsCellsN3 =  deepcopy( calculateVariablesAtStage(1.0/2.0,UconsCellsN2,UphysCells, UphysNodes) );
		# #UconsCellsNew = deepcopy( calculateVariablesAtStage(1.0,UconsCellsOld,UphysCells, UphysNodes) );
		
		
		# UconsCellsNew = deepcopy(FirstOrderUpwindM1(1.0,UconsCellsOld, UphysCells ));
		# #UconsCellsNew = deepcopy( SecondOrderUpwind(1.0,UconsCellsOld, UphysCells, UphysNodes));
		
	# else

		# display("parallel solution is not implemented");
		
		# ################################################################
		# ### Euler first order
		# #calculateVariablesAtStageX11(1.0,UconsCellsOldX,UconsCellsNewX);
		# #SecondOrderUpwindParallelX11(1.0,UconsCellsOldX,UconsCellsNewX);
		
		# ################################################################
		# ### RK4 
		# #calculateVariablesAtStageX11(1.0/4.0,UconsCellsOldX,UconsCellsN1X);
		# #calculateVariablesAtStageX11(1.0/3.0,UconsCellsN1X,UconsCellsN2X);
		# #calculateVariablesAtStageX11(1.0/2.0,UconsCellsN2X,UconsCellsN3X);
		# #calculateVariablesAtStageX11(1.0,UconsCellsN3X,UconsCellsNewX);
		
	# end

# elseif (solver.TimeDiscretization == 0) 

	# #const rk21  = 1.0;
	# #const rk22  = 0.5;

	# #UconsCellsN1 = calculateVariablesAtStage(rk21,UconsCellsOld,thermo.RGAS);

	# #UconsCellsN2 = UconsCellsOld*rk22 + UconsCellsN1*rk22;
	# #UconsCellsNew = calculateVariablesAtStage(rk22,UconsCellsN2,thermo.RGAS); 	
# end









