

include("RoeFlux2d.jl")
include("AUSMflux2d.jl"); #AUSM+ inviscid flux calculation 
include("utilsFVM2.jl"); #FVM utililities
include("RIEMANN1d.jl")

include("thermo.jl"); #setup thermodynamics



phsLeft = zeros(Float64,4);
phsRight = zeros(Float64,4);

phsLeft[1] = 1.0;
phsLeft[2] = 290.0;
phsLeft[3] = 0.0;
phsLeft[4] = 7143.0;

phsRight[1] = 1.7;
phsRight[2] = 263.72;
phsRight[3] = -51.62;
phsRight[4] = 15282.0;


nx = 1.0;
ny = 0.0;

side = 1.0;


fluxTmpRoe = zeros(Float64,4);
fluxTmpAUSM = zeros(Float64,4);
fluxTmpRiemann= zeros(Float64,4);

fluxTmpRoe = RoeFlux2d(phsRight,phsLeft, nx,ny,side,thermo.Gamma);
fluxTmpAUSM = AUSMplusFlux2d(phsRight,phsLeft, nx,ny,side,thermo.Gamma);

fluxTmpRiemann = compute_1D_ARBITRARY_INVISCID_RIEMANN_FLUX_from_UPHYS(phsRight,phsLeft,nx,ny,side,thermo.Gamma);

println("Calculated fluxes: \n");
println("Roe: \t", fluxTmpRoe[1], "\tAUSM: \t",fluxTmpAUSM[1],"\tRiemann: \t", fluxTmpRiemann[1] );
println("Roe: \t", fluxTmpRoe[2], "\tAUSM: \t",fluxTmpAUSM[2],"\tRiemann: \t", fluxTmpRiemann[2] );
println("Roe: \t", fluxTmpRoe[3], "\tAUSM: \t",fluxTmpAUSM[3],"\tRiemann: \t", fluxTmpRiemann[3] );
println("Roe: \t", fluxTmpRoe[4], "\tAUSM: \t",fluxTmpAUSM[4],"\tRiemann: \t", fluxTmpRiemann[4] );