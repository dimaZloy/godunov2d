Solve 2d Euler equations based on the Godunov finite-volume method (FVM).
Available space approximation schemes are: 
the first-order (FOU) and second-order upwind with the Minmod flux limiter (SOU).
The inviscid flux is calculated using the AUSMplus and Roe methods. 
Time marching is done using the explicit Euler method. 
Utilize Paraview for post-processing. 

Tested on a steady-state shock wave reflecting from a rigid surface, published originally   
by Yee et al (1984), "Implicit Total Variation Diminishing (TVD) Schemes for Steady-State Calculations," AIAA paper 83-1902. 

The computaion domain is a rectangle [0 4.1]x [0 1] discretized by 80x20 cells. 
The following boundary conditions are applied (assuming perfect gas and Gamma = 1.4):
Left: rho = 1.0; P = 7143; U = 290; V = 0.0; 
Top: rho = 1.7; P = 15282; U = 261.193; V = -0.50632;
Bottom: solid wall;
Right: zero-gradient (supersonic flow). 
Simple triangular, quad and mixed grids are available for testing.   

The analytical solution for this oblique shock problem (with ksi = 29 deg) 
leeds to the reflection angle of 23.279 deg and flow parameters (at the right boundary): 
rho = 2.6872; P = 29340; U = 2.4015;  V = 0. 

![alt tag](https://github.com/dimaZloy/godunov2d/blob/main/results1.png) 
