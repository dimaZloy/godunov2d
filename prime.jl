

using PyPlot;
using WriteVTK;
using CPUTime;
using Distributed;
using SharedArrays;
using DelimitedFiles;
using Printf
using BSON: @load
using BSON: @save


include("D:/Julia/JuliaProjects/aero2d/mesh2d/primeObjects.jl");
include("thermo.jl"); #setup thermodynamics
include("utilsIO.jl");
include("limiters.jl");
include("RoeFlux2d.jl")
include("AUSMflux2d.jl"); #AUSM+ inviscid flux calculation 
include("utilsFVM2d.jl"); #FVM utililities

include("computeslope2d.jl")
include("SOUscheme.jl")
include("createFields2d.jl");


## TODO: 
## check:  how to usemaxArea or maxEdgeLength

##TODO:
## test tri-mesh2d
## test quad mesh2d
## test mixed tri-quad



function godunov2d(pname::String, approx::String, outputfile::String)

	
	#@load "testTriMesh2d.bson" testMesh
	#@load "testQuadMesh2d.bson" testMesh
	#@load "testMixedMesh2d.bson" testMesh
	
	@load pname testMesh

	include("setupSolver2DoblickShock2d.jl"); #setup FVM and numerical schemes
	
		
	# display(testMesh.node2cellsL2up)
	# display(testMesh.node2cellsL2down)

	
	timeVector = [];
	residualsVector1 = []; 
	residualsVector2 = []; 
	residualsVector3 = []; 
	residualsVector4 = []; 

	#Delta2 = zeros(Float64,testMesh.nCells,4); 
	residualsVectorMax = ones(Float64,4);
	convergenceCriteria= [1e-5;1e-5;1e-5;1e-5;];

	## init primitive variables 
	println("set initial and boundary conditions ...");
	testfields2d = createFields2d(testMesh, thermo);
	
	# debug = true;
	# if (debug)
	
		## init conservative variables 
		UconsCellsOld = zeros(Float64,testMesh.nCells,4);
		UconsCellsNew = zeros(Float64,testMesh.nCells,4);
		Delta = zeros(Float64,testMesh.nCells,4);

		UconsCellsOld = phs2dcns2dcells(testfields2d, thermo.Gamma); #old vector
		UconsCellsNew = phs2dcns2dcells(testfields2d, thermo.Gamma); #new vector
		
		
		# UconsCellsN1 = deepcopy(UconsCellsOld); #new  vector 
		# UconsCellsN2 = deepcopy(UconsCellsOld); #new  vector 
		# UconsCellsN3 = deepcopy(UconsCellsOld); #new  vector 
		
		
		
		(zzz2,id) = findmin(testMesh.cell_edges_length[:,1]);
		(zzz3,id) = findmin(testMesh.cell_edges_length[:,2]);
		(zzz4,id) = findmin(testMesh.cell_edges_length[:,3]);
		(zzz5,id) = findmin(testMesh.cell_edges_length[:,4]);
		
		DX::Float64 = 0.25*sqrt(zzz2*zzz2 + zzz3*zzz3 + zzz4*zzz4 + zzz5*zzz5);
				
		
		println("Start calculations ...");
		println(output.header);

		while (dynControls.isRunSimulation == 1)
			CPUtic();	
			
			
			# PROPAGATE STAGE: 
			(dynControls.velmax,id) = findmax(testfields2d.VMAXCells);
			#dynControls.tau = 0.5*solControls.CFL * testMesh.maxArea/(max(dynControls.velmax,1.0e-6));
			dynControls.tau =  solControls.CFL * DX/(max(dynControls.velmax,1.0e-6));
			
						
			dt::Float64 = 0.0;
			if (solControls.timeStepMethod == 1)
				dt  = dynControls.tau;  
			else
				dt =  solControls.dt;  
			end
			
			#SecondOrderUpwindM2(1.0, UconsCellsOld, node2cellL2up, node2cellL2down, testMesh, testfields2d, thermo, solControls, dynControls, UconsCellsNew );
			SecondOrderUpwindM2(1.0, dt,  UconsCellsOld, testMesh, testfields2d, thermo, UconsCellsNew );
			
			# rk21  = 1.0;
			# rk22  = 0.5;

			# UconsCellsN1 = deepcopy(SecondOrderUpwindM2(rk21, UconsCellsOld, node2cellL2up, node2cellL2down, testMesh, testfields2d, thermo, solControls, dynControls ));
			# UconsCellsN2 = deepcopy(UconsCellsOld.*rk22 .+ UconsCellsN1.*rk22);
			# UconsCellsNew = deepcopy(SecondOrderUpwindM2(rk22, UconsCellsN2, node2cellL2up, node2cellL2down, testMesh, testfields2d, thermo, solControls, dynControls )); 	
		
		
		
			# EVALUATE STAGE:
			
			dynControls.flowTime += dt; 
			
			push!(timeVector, dynControls.flowTime); 
			dynControls.curIter += 1; 
			dynControls.verIter += 1;
			
			updateVariablesM2(Delta, UconsCellsOld, UconsCellsNew,  testMesh, testfields2d, dynControls);
			#updateVariables!(UconsCellsOld, testMesh, testfields2d, dynControls);
						
			
			updateResidual!(Delta, 
				residualsVector1,residualsVector2,residualsVector3,residualsVector4, residualsVectorMax,  
				convergenceCriteria, dynControls);
			
			updateOutput!(timeVector,residualsVector1,residualsVector2,residualsVector3,residualsVector4, residualsVectorMax, 
				testMesh, testfields2d, solControls, output, dynControls);
			

			if (flowTime>= solControls.stopTime || dynControls.isSolutionConverged == 1)
				dynControls.isRunSimulation = 0;
		
				if (dynControls.isSolutionConverged == true)
					println("Solution converged! ");
				else
					println("Simultaion flow time reached the set Time!");
				end
			
				if (output.saveResiduals == 1)
					println("Saving Residuals ... ");
					cd(dynControls.localTestPath);
					saveResiduals(output.fileNameResiduals, timeVector, residualsVector1, residualsVector2, residualsVector3, residualsVector4);
					cd(dynControls.globalPath);
				end
				if (output.saveResults == 1)
					println("Saving Results ... ");
					cd(dynControls.localTestPath);
					saveSolution(output.fileNameResults, testMesh.xNodes, testMesh.yNodes, UphysNodes);
					cd(dynControls.globalPath);
				end
			
			end

			dynControls.cpuTime  += CPUtoq(); 
			
			if (dynControls.flowTime >= solControls.stopTime)
				dynControls.isRunSimulation = 0;
			end
			
		end ## end while
		
		
		
		println("Saving  solution to  ", outputfile);
		saveResults2VTK(outputfile, testMesh, testfields2d.densityNodes, "density");
		println("done ...  ");	
		 
	#end ## end debug
	
end


# @time godunov2d("testTriMesh2d.bson","first","fouTri"); 
# @time godunov2d("testQuadMesh2d.bson","first","fouQuad"); 
# @time godunov2d("testMixedMesh2d.bson","first","fouMixed"); 

@time godunov2d("testTriMesh2d.bson","second","souTri"); 
#@time godunov2d("testQuadMesh2d.bson","second","souQuad"); 
#@time godunov2d("testMixedMesh2d.bson","second","souMixed"); 





