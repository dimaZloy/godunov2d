

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
include("RoeFlux2d.jl")
include("AUSMflux2d.jl"); #AUSM+ inviscid flux calculation 
include("utilsFVM2d.jl"); #FVM utililities
include("RIEMANN1d.jl")


include("createFields2d.jl");

## TODO: 
## check:  how to usemaxArea or maxEdgeLength
## check:  cells2nodesSolutionReconstructionWithStencils(testMesh, testFields.densityCells ) with testfields2d structure 

function godunov2d()

	
	@load "testMesh02.bson" testMesh
	
	testMesh;
	
	
	#display(testMesh.maxArea)
	#display(sqrt(testMesh.maxEdgeLength))

	
	global testdir = "D:/Julia/JuliaProjects/aero2d/godunov2d";

	include("setupSolver2DoblickShock2d.jl"); #setup FVM and numerical schemes
	
	
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
	
	
	## init conservative variables 
	UconsCellsOld = zeros(Float64,testMesh.nCells,4);
	UconsCellsNew = zeros(Float64,testMesh.nCells,4);
	Delta = zeros(Float64,testMesh.nCells,4);

	UconsCellsOld = phs2dcns2dcells(testfields2d, thermo.Gamma); #old vector
	UconsCellsNew = deepcopy(UconsCellsOld); #new  vector 
	
	println("Start calculations ...");
	println(output.header);

	 while (dynControls.isRunSimulation == 1)
		CPUtic();	
		
		
		# PROPAGATE STAGE: 
		(dynControls.velmax,id) = findmax(testfields2d.VMAXCells);
		#dynControls.tau = solControls.CFL * testMesh.maxEdgeLength/(max(dynControls.velmax,1.0e-6));
		dynControls.tau = solControls.CFL * testMesh.maxArea/(max(dynControls.velmax,1.0e-6));
		
		
		#UconsCellsNew = deepcopy(FirstOrderUpwindM2(1.0, UconsCellsOld, testfields2d ));
		UconsCellsNew = deepcopy(FirstOrderUpwindM2(1.0, UconsCellsOld, testMesh, testfields2d, thermo, solControls, dynControls ));
		

	
		# EVALUATE STAGE:
		if (solControls.timeStepMethod == 1)
			dynControls.flowTime += dynControls.tau;  
		else
			dynControls.flowTime += solControls.dt;  
		end
		
		push!(timeVector, dynControls.flowTime); 
		dynControls.curIter += 1; 
		dynControls.verIter += 1;

		
		Delta =  deepcopy(UconsCellsNew - UconsCellsOld); 
		UconsCellsOld =  deepcopy(UconsCellsNew);	
		
		updateVariables!(UconsCellsOld, testMesh, testfields2d, dynControls);
		
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
		
	 end
	
end


@time godunov2d(); 





