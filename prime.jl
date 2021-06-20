

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
include("utilsIO.jl");
include("RoeFlux2d.jl")
include("AUSMflux2d.jl"); #AUSM+ inviscid flux calculation 
include("utilsFVM2d.jl"); #FVM utililities
include("RIEMANN1d.jl")

include("thermo.jl"); #setup thermodynamics

## TODO: 
## make mutable structure for uphys and ucons 


function godunov2d()

	
	@load "testMesh02.bson" testMesh
	
	global testMesh;

	## create VTK cell types to save solution to VTK
	## VTKcells = createVTKcells(testMesh); 
	
	global testdir = "D:/Julia/JuliaProjects/aero2d/godunov2d";

	include("setupSolver2DoblickShock2d.jl"); #setup FVM and numerical schemes
	include("setupTest2DoblickShock2d.jl"); #setup boundary and initial conditions


	while (dynControls.isRunSimulation == 1)
		CPUtic();	
		include("propagate.jl");
		include("evaluate.jl");
	
		dynControls.cpuTime  += CPUtoq(); 
	end
	
end


godunov2d(); 





