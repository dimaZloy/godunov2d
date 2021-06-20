@inline  @everywhere  function pow(x::Float64, y::Float64)::Float64
	return x^y; 
end


@inline  @everywhere  function dfgod(pBig::Float64, p::Float64,rho::Float64, gamma::Float64)::Float64

	Tmp1_z::Float64 = (gamma+1.0)/(2.0*gamma);
	Tmp2_z::Float64 = (gamma-1.0)/(2.0*gamma);
	Tmp3_z::Float64 = 2.0/(gamma-1.0);
	Tmp4_z::Float64 = (gamma-1.0)/2.0/gamma;

	pRel::Float64 = pBig/p;
	cSound::Float64=sqrt(gamma*p/rho);
	dp::Float64=pBig-p;

	f1::Float64=((gamma+1.0)*pRel+(3.0*gamma-1.0))/(4.0*gamma*rho*cSound*pow(Tmp1_z*pRel+Tmp2_z,1.5));
	f2::Float64=(1.0/gamma*pBig)*cSound*pow(pRel,Tmp2_z);   
			
	return (f1-f2)*(dp>0)+f2 ;
end

@inline  @everywhere  function fgod( pBig::Float64, p::Float64,rho::Float64, gamma::Float64)::Float64

	Tmp1_z::Float64 = (gamma+1.0)/(2.0*gamma);
	Tmp2_z::Float64 = (gamma-1.0)/(2.0*gamma);
	Tmp3_z::Float64 = 2.0/(gamma-1.0);
	Tmp4_z::Float64 = (gamma-1.0)/2.0/gamma;	

	pRel::Float64 = pBig/p;
	cSound::Float64 = sqrt(gamma*p/rho);
	dp::Float64 = pBig-p;
	
	f1::Float64 = dp/( rho*cSound* sqrt(Tmp1_z*pRel+Tmp2_z) );
	f2::Float64 = Tmp3_z*cSound*(pow(pRel,Tmp2_z)-1.0);
	
	return (f1-f2)*(dp>0) + f2;
end


@inline  @everywhere  function computeRiemannFlux(rhoLeft::Float64,uLeft::Float64,pLeft::Float64,
		rhoRight::Float64, uRight::Float64,pRight::Float64,  Gamma::Float64)::Array{Float64,1}


	RS_EpsMethod = 1.0e-6;
	RS_IterationLimit = 25;
	kgas::Float64 = (Gamma-1.0)/(Gamma+1.0);
	

	cSoundLeft::Float64  = sqrt(Gamma*pLeft/rhoLeft);
	cSoundRight::Float64 = sqrt(Gamma*pRight/rhoRight);

	
	if(abs(pLeft-pRight)<0.01)
		pRight = pRight +  1.0;
	end
	
	
	pMin::Float64  = min(pLeft,pRight);
	pMax::Float64  = max(pLeft,pRight);
	zero::Float64 = 0.0;


	VShock::Float64  = fgod(pMax,pLeft,rhoLeft,Gamma) + fgod(pMax,pRight,rhoRight,Gamma);
	VRar::Float64  = fgod(pMin,pLeft,rhoLeft,Gamma) + fgod(pMin,pRight,rhoRight,Gamma);
	VVac::Float64  = fgod(zero,pLeft,rhoLeft,Gamma)+ fgod(zero,pRight,rhoRight,Gamma);
	
	Degree::Float64  = log((VShock-VVac)/(VRar-VVac))/log(pMax/pMin);
	Coeff::Float64   = (VShock-VVac)/pow(pMax,Degree);



	#pBig::Float64  = pow( ((uLeft-uRight-VVac)/Coeff), (1.0/Degree) );

	pBig::Float64  = ((uLeft-uRight-VVac)/Coeff)^(1.0/Degree);
   
	niter::Int64 = 0;


	while (true)
	
		
		niter = niter + 1;
		fg1::Float64 = fgod(pBig,pLeft,rhoLeft,Gamma);
		fg2::Float64 = fgod(pBig,pRight,rhoRight,Gamma);
		dfg1::Float64 = dfgod(pBig,pLeft,rhoLeft,Gamma);
		dfg2::Float64 = dfgod(pBig,pRight,rhoRight,Gamma);
		delta::Float64 =  (fg1 + fg2 - (uLeft-uRight))/( dfg1 + dfg2);
		pBigNew::Float64 = pBig - delta; 
		ThisEpsilon::Float64 =  abs(pBigNew-pBig);
		pBig = pBigNew;
		
		if ((abs(ThisEpsilon)<=RS_EpsMethod)||(niter>=RS_IterationLimit))
			break;
		end
		##if ((abs(ThisEpsilon)<=RS_EpsMethod)&&(niter>=RS_IterationLimit))
		##	break;
		##end
		
	end
	
	
	dV::Float64 = uLeft-uRight;

	CondP1b = dV>VShock;
	CondP2b = (VRar<dV)&(dV<=VShock)&(pRight>pLeft); 
	CondP3b = (VRar<dV)&(dV<=VShock)&(pRight<=pLeft);
	CondP4b = (VVac<dV)&(dV<=VRar);
	CondP5b = dV<=VVac;



	aLeftShock::Float64 = sqrt(rhoLeft*((Gamma+1.0)/2.0*pBig + (Gamma-1.0)/2.0*pLeft));
                       
	aLeftFan::Float64 = ((Gamma-1.0)/(2.0*Gamma))*rhoLeft*cSoundLeft*(1.0-(pBig)/(pLeft))/(1.0-pow(abs(pBig/pLeft),(Gamma-1)/(2*Gamma)) + eps() );

	if (abs(pBig-pLeft)<0.01)
		   aLeftFan = rhoLeft*cSoundLeft;
	end
	
	tmpb  = CondP1b|CondP2b;
	tmpb1 = CondP3b|CondP4b|CondP5b;


	aLeft::Float64 = aLeftShock*tmpb + aLeftFan*tmpb1;
    
	aRightShock::Float64 = sqrt(rhoRight*((Gamma+1.0)/2.0*pBig + (Gamma -1.0)/2.0*pRight));

        aRightFan::Float64 = ((Gamma-1.0)/(2.0*Gamma))*rhoRight*cSoundRight*(1.0-(pBig)/(pRight))/(1.0- pow(abs(pBig/pRight),((Gamma-1.0)/(2.0*Gamma))) + eps() );   


	if (abs(pBig-pRight)<0.01)
		aRightFan = rhoRight*cSoundRight;
	end


	tmpb  = CondP1b|CondP3b;
	tmpb1 = CondP2b|CondP4b|CondP5b;


    	aRight::Float64 = aRightShock*tmpb + aRightFan*tmpb1;


	uBig::Float64 = (aLeft*uLeft + aRight*uRight+pLeft-pRight)/(aLeft+aRight)*!CondP5b;

    	
	rhoBigLeftShock::Float64 = rhoLeft*aLeft/(aLeft-rhoLeft*(uLeft-uBig));
	cSoundLeftPrim::Float64 = cSoundLeft + 0.5*(Gamma-1.0)*(uLeft-uBig);
	rhoBigLeftFan::Float64 = Gamma*pBig/pow(cSoundLeftPrim,2.0);


	tmpb  = CondP1b|CondP2b;
	tmpb1 = CondP3b|CondP4b;

	rhoBigLeft::Float64 = rhoBigLeftShock*tmpb +  rhoBigLeftFan*tmpb1;

	
	rhoBigRightShock::Float64 = rhoRight*aRight/(aRight+rhoRight*(uRight-uBig));
	cSoundRightPrim::Float64 = cSoundRight-0.5*(Gamma-1.0)*(uRight-uBig);
	rhoBigRightFan::Float64 = Gamma*pBig/pow(cSoundRightPrim,2.0);


	tmpb  = CondP1b|CondP3b;
	tmpb1 = CondP2b|CondP4b;

	rhoBigRight::Float64 = rhoBigRightShock*tmpb + rhoBigRightFan*tmpb1;
    
	
	uLeftShock::Float64 = uLeft - aLeft/rhoLeft;
	uLeftFan::Float64 = uLeft - cSoundLeft;
	uLeftFanBack::Float64 = uBig- cSoundLeftPrim;


	tmpb  = CondP1b|CondP2b;
	tmpb1 = CondP3b|CondP4b|CondP5b;

	uLeftFront::Float64 = uLeftShock*tmpb + uLeftFan*tmpb1;
	uLeftBack::Float64  = uLeftShock*tmpb + uLeftFanBack*tmpb1;
	
    
	uRightShock::Float64   = uRight +  aRight/rhoRight;
	uRightFan::Float64     = uRight + cSoundRight;
	uRightFanBack::Float64 = uBig + cSoundRightPrim;

	tmpb  = CondP1b|CondP3b;
	tmpb1 = CondP2b|CondP4b|CondP5b;


	uRightFront::Float64 = uRightShock*tmpb + uRightFan*tmpb1;
	uRightBack::Float64  = uRightShock*tmpb + uRightFanBack*tmpb1;



	tmpb1  = uBig>=0.0;
	tmpb2 = uLeftFront>0.0;
	tmpb3 = uLeftBack<0.0;
	tmpb4 = (uLeftFront<0.0)&(uLeftBack>=0.0);
	tmpb5 = uBig<0.0;
	tmpb6 = uRightFront<0.0;
	tmpb7 = uRightBack>0.0;
	tmpb8 = (uRightFront>0.0)&(uRightBack<=0.0);


	uSide::Float64=tmpb1*(tmpb2*uLeft + tmpb3*uBig + tmpb4*(kgas*(2.0/(Gamma-1.0)*cSoundLeft+uLeft)))+ tmpb5*(tmpb6*uRight + tmpb7*uBig + tmpb8*(kgas*(uRight-2.0/(Gamma-1.0)*cSoundRight)));


	rhoSide::Float64 = tmpb1*(tmpb2*rhoLeft + tmpb3*rhoBigLeft + tmpb4*rhoLeft*pow(uSide/cSoundLeft,5.0)) +	 tmpb5*(tmpb6*rhoRight + tmpb7*rhoBigRight + tmpb8*rhoRight*pow(-uSide/cSoundLeft,5.0));

 	
	pSide::Float64 = tmpb1*(tmpb2*pLeft + tmpb3*pBig + tmpb4*rhoSide * uSide * uSide/Gamma)  +tmpb5*(tmpb6*pRight + tmpb7*pBig +  tmpb8*rhoSide* uSide * uSide/Gamma);
    
	

	speeds = zeros(Float64,6); 
	speeds[1] = abs(uLeftShock);
	speeds[2] = abs(uLeftFan);
	speeds[3] = abs(uBig-cSoundLeftPrim);
	speeds[4] = abs(uRightShock);
	speeds[5] = abs(uRightFan);
	speeds[6] = abs(uBig+cSoundRightPrim);

	maxV,dummy1  = findmax(speeds);          

	vSide = zeros(Float64,4) 

	#vSide[1] = rhoSide*uSide;
	#vSide[2] = rhoSide*uSide*uSide + pSide;
	#vSide[3] = ( Gamma/(Gamma-1.0)*pSide + 0.5*rhoSide*uSide*uSide)*uSide;
	

	#Flux=[rhoSide.*uSide; ...
	#      rhoSide.*uSide.^2+pSide; ...
        #(Gamma/(Gamma-1)*pSide+0.5*rhoSide.*uSide.^2).*uSide];
	
	vSide[1] = rhoSide;
	vSide[2] = uSide;
	vSide[3] = pSide;


	
	vSide[4] = maxV;


	return vSide; 

end


# @inline  function compute_1D_ARBITRARY_INVISCID_RIEMANN_FLUX_from_UPHYS(uLeft::Array{Float64,2}, uRight::Array{Float64,2}, nx::Float64,  ny::Float64, side::Float64, gamma::Float64)
	# return compute_1D_ARBITRARY_INVISCID_RIEMANN_FLUX_from_UPHYS(uLeft[1],uLeft[2],uLeft[3],uLeft[4], uRight[1], uRight[2],uRight[3],uRight[4], nx,ny,side, gamma); 
# end

@inline @everywhere  function RiemannFlux2d(
			uLeft::Array{Float64,1}, uRight::Array{Float64,1}, 
			nx::Float64,  ny::Float64, side::Float64, gamma::Float64)::Array{Float64,1}
			
	return compute_1D_ARBITRARY_INVISCID_RIEMANN_FLUX_from_UPHYS(
		uLeft[1],uLeft[2],uLeft[3],uLeft[4], uRight[1], uRight[2],uRight[3],uRight[4], nx, ny,side, gamma); 
end


@inline @everywhere   function compute_1D_ARBITRARY_INVISCID_RIEMANN_FLUX_from_UPHYS(
				rhoL::Float64,	UL::Float64, 	VL::Float64, 	PL::Float64,
				rhoR::Float64,	UR::Float64, 	VR::Float64,	PR::Float64, 
				nx::Float64,  ny::Float64, side::Float64, 
				gamma::Float64)::Array{Float64,1}

		flux = zeros(Float64,4);


		NLeft::Float64  = UL*nx + VL*ny;
		TLeft::Float64  = UL*ny - VL*nx;

		NRight::Float64 = UR*nx + VR*ny;
		TRight::Float64 = UR*ny - VR*nx;

	
		zvar = zeros(Float64,4);
		
		zvar = computeRiemannFlux(rhoL,NLeft,PL,rhoR,NRight,PR,gamma);

		R::Float64  = zvar[1];
		Vn::Float64 = zvar[2];
		P::Float64 = zvar[3];
		Vmax::Float64 = zvar[4];
		
		TLeft = 0.5*(TLeft+TRight);
		#TLeft = (TLeft+TRight);


		U::Float64 = Vn*nx + TLeft*ny;
		V::Float64 = Vn*ny - TLeft*nx;

		
		E::Float64 = P/(gamma-1.0)+0.5*R*(U*U+V*V);
		N::Float64 = U*nx + V*ny;
 
		# display(TLeft)
		# display(U)
		# display(V)
		# display(E)
		# display(N)
 
		flux[1] = -(R*N) * side;
		flux[2] = -(R*N*U + nx*P) * side;
		flux[3] = -(R*N*V + ny*P) * side;
		flux[4] = -(R*N*E ) * side;
		
		return flux;

end



@inline @everywhere  function GOD1d(NCells::Int64, UPHSCELLS::Array{Float64,2}, Gamma::Float64)::Array{Float64,2}

	UFLUX = zeros(NCells,4);
	rhoL::Float64 = 0.0;
	uL::Float64 = 0.0;
	pL::Float64 = 0.0;

	rhoR::Float64 = 0.0;
	uR::Float64 = 0.0;
	pR::Float64 = 0.0;

	flux = zeros(Float64,4);

	@simd for i=1:NCells-1
	   rhoL = UPHSCELLS[i,1];
	   uL = UPHSCELLS[i,2];
   	   pL = UPHSCELLS[i,3];	

	   rhoR = UPHSCELLS[i+1,1];
	   uR = UPHSCELLS[i+1,2];
   	   pR = UPHSCELLS[i+1,3];	

	   flux = computeRiemannFlux(rhoL,  uL, pL, rhoR, uR, pR, Gamma);
	  ##  flux = computeRiemannFlux( rhoR, uR, pR, rhoL,  uL, pL, Gamma);
	   ##flux = computeRiemannFlux2(rhoL,  uL, 0.0, pL, rhoR, uR, 0.0, pR, 1.0,0.0, 1.0, Gamma);
	  
           UFLUX[i,1:4] = flux;
				
	end

	return UFLUX;
	
end





