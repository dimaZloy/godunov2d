
println("set thermo-physical properties ...");

@everywhere struct THERMOPHYSICS
	RGAS::Float64;
	Gamma::Float64;
	Cp::Float64;
	mu::Float64;
	kGas::Float64;		
end

thermo = THERMOPHYSICS(287.0,1.4,1000.0,0.0,0.0);

#@everywhere const thermoX = $thermo;




