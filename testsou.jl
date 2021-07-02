
# @everywhere struct mesh2d
	# nCells::Int64
	# nNodes::Int64
	# nNeibCells::Int64						## max number of neighbors 
	# nBSets::Int64							##  number of boundaries  
	# xNodes::Array{Float64,1} 				##  mesh_nodes[nNodesx3]
	# yNodes::Array{Float64,1} 				##	mesh_nodes[nNodesx3]
	# mesh_connectivity::Array{Int64,2} 	## [nCellsx3]
	# bc_data::Array{Int64,2}
	# bc_indexes::Array{Int64,1}
	# cell_nodes_X::Array{Float64,2} 			## [nCellsx4]
	# cell_nodes_Y::Array{Float64,2} 			## [nCellsx4]
	# cell_mid_points::Array{Float64,2} 		## [nCellsx2]
	# cell_areas::Array{Float64,1} 			## [nCellsx1]
	# Z::Array{Float64,1} 					## [nCellsx1] 1/cell_areas
	# cell_edges_Nx::Array{Float64,2} 		## [nCellsx4]
	# cell_edges_Ny::Array{Float64,2} 		## [nCellsx4]
	# cell_edges_length::Array{Float64,2} 	## [nCellsx4]
	# cell_stiffness::Array{Int64,2} 		## [nCellsx4]
	# cell_clusters::Array{Int64,2} 		## [nNodesx8]
	# node_stencils::Array{Float64,2} 		## [nNodesx8]
	# maxArea::Float64
	# maxEdgeLength::Float64
	# VTKCells::Array{MeshCell,1}
# end

# mesh_connectivity has the following format (nCells x 7):
# 1 element - global id,
# 2 element - cell type (2 for quads /3 for triangles)
# 3 element - number of nodes
# 4-7 - nodes ids
# 3 nodes for tri mesh 
# 4 nodes for quad mesh 


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
include("computeNode2CellsL2.jl")



@load "testTriMesh2d.bson" testMesh

	
	
node2cellL2up = zeros(Int64,testMesh.nCells,4);
node2cellL2down = zeros(Int64,testMesh.nCells,4);
	
computeNode2CellsL2(testMesh, node2cellL2up, node2cellL2down);


N = 901;

el_type =   testMesh.mesh_connectivity[N,2];
num_nodes = testMesh.mesh_connectivity[N,3];

elementXleft = [];
elementYleft = [];

push!(elementXleft,testMesh.cell_nodes_X[N,1])
push!(elementXleft,testMesh.cell_nodes_X[N,2])
push!(elementXleft,testMesh.cell_nodes_X[N,3])
push!(elementXleft,testMesh.cell_nodes_X[N,1])


push!(elementYleft,testMesh.cell_nodes_Y[N,1])
push!(elementYleft,testMesh.cell_nodes_Y[N,2])
push!(elementYleft,testMesh.cell_nodes_Y[N,3])
push!(elementYleft,testMesh.cell_nodes_Y[N,1])



neib_cell_1 = testMesh.cell_stiffness[N,1];
neib_cell_2 = testMesh.cell_stiffness[N,2];
neib_cell_3 = testMesh.cell_stiffness[N,3];

elRightX1 = [];
elRightX2 = [];
elRightX3 = [];

elRightY1 = [];
elRightY2 = [];
elRightY3 = [];



push!(elRightX1,testMesh.cell_nodes_X[neib_cell_1,1])
push!(elRightX1,testMesh.cell_nodes_X[neib_cell_1,2])
push!(elRightX1,testMesh.cell_nodes_X[neib_cell_1,3])
push!(elRightX1,testMesh.cell_nodes_X[neib_cell_1,1])

push!(elRightY1,testMesh.cell_nodes_Y[neib_cell_1,1])
push!(elRightY1,testMesh.cell_nodes_Y[neib_cell_1,2])
push!(elRightY1,testMesh.cell_nodes_Y[neib_cell_1,3])
push!(elRightY1,testMesh.cell_nodes_Y[neib_cell_1,1])


push!(elRightX2,testMesh.cell_nodes_X[neib_cell_2,1])
push!(elRightX2,testMesh.cell_nodes_X[neib_cell_2,2])
push!(elRightX2,testMesh.cell_nodes_X[neib_cell_2,3])
push!(elRightX2,testMesh.cell_nodes_X[neib_cell_2,1])

push!(elRightY2,testMesh.cell_nodes_Y[neib_cell_2,1])
push!(elRightY2,testMesh.cell_nodes_Y[neib_cell_2,2])
push!(elRightY2,testMesh.cell_nodes_Y[neib_cell_2,3])
push!(elRightY2,testMesh.cell_nodes_Y[neib_cell_2,1])


push!(elRightX3,testMesh.cell_nodes_X[neib_cell_3,1])
push!(elRightX3,testMesh.cell_nodes_X[neib_cell_3,2])
push!(elRightX3,testMesh.cell_nodes_X[neib_cell_3,3])
push!(elRightX3,testMesh.cell_nodes_X[neib_cell_3,1])


push!(elRightY3,testMesh.cell_nodes_Y[neib_cell_3,1])
push!(elRightY3,testMesh.cell_nodes_Y[neib_cell_3,2])
push!(elRightY3,testMesh.cell_nodes_Y[neib_cell_3,3])
push!(elRightY3,testMesh.cell_nodes_Y[neib_cell_3,1])




upnode1 = node2cellL2up[N,1];
downnode1 = node2cellL2down[N,1];

upnode2 = node2cellL2up[N,2];
downnode2 = node2cellL2down[N,2];

upnode3 = node2cellL2up[N,3];
downnode3 = node2cellL2down[N,3];



figure(1)
clf()

plot(elementXleft, elementYleft, "-k",linewidth = 1.0);

plot(elRightX1, elRightY1, "--r",linewidth = 1.0);
plot(elRightX2, elRightY2, "--g",linewidth = 1.0);
plot(elRightX3, elRightY3, "--b",linewidth = 1.0);



plot(testMesh.xNodes[upnode1],testMesh.yNodes[upnode1],"or",markersize = 5.0);
plot(testMesh.xNodes[downnode1],testMesh.yNodes[downnode1],"sr",markersize = 7.0);

plot(testMesh.xNodes[upnode2],testMesh.yNodes[upnode2],"og",markersize = 5.0);
plot(testMesh.xNodes[downnode2],testMesh.yNodes[downnode2],"sg",markersize = 7.0);

plot(testMesh.xNodes[upnode3],testMesh.yNodes[upnode3],"ob",markersize = 5.0);
plot(testMesh.xNodes[downnode3],testMesh.yNodes[downnode3],"sb",markersize = 7.0);


#tricontourf(testMesh.xNodes,testMesh.yNodes, testFields.densityNodes,pControls.nContours,vmin=pControls.rhoMINcont,vmax=pControls.rhoMAXcont);
#tricontour(testMesh.xNodes,testMesh.yNodes, testFields.densityNodes,pControls.nContours,vmin=pControls.rhoMINcont,vmax=pControls.rhoMAXcont);
			
			
#xlim(0.0,4.1)
#ylim(0.0,1.0)

xlabel("x");
ylabel("y");
#title("Contours of density");
#axis("equal");