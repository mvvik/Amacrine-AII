function ModelOutput = GetModelOutput5( EQ, CaGrid, timeArray, rArray, X, DT, odeSolve, ACC )

%   Parameters: distances, sensor binding kinetics, and pool fraction:
%      1      2      3        4       5    6        7   
%   Dist-1 Dist-2 SensorKD SensorKp coop gamma pool-fraction
    
    totalTime   = DT(end);
	[TT, RR]    = meshgrid(timeArray, rArray);
			
    KD          = X(3);  
    kp          = X(4);  
    km          = KD * kp;
	coop        = X(5);
	gamma       = X(6);

	kp1 = 5*kp;        kp2 = 4*kp;      kp3 = 3*kp; kp4 = 2*kp;      kp5 = kp;    
    km1 = km/(coop^2); km2 = 2*km/coop; km3 = 3*km; km4 = 4*km*coop; km5 = 5*km*coop^2;

    fraction    = [ X(7),  1 - X(7) ];
	ModelOutput = zeros( size( DT ) );

    for Iter = [1 2]
	
		R        = X( Iter );
		[tt, rr] = meshgrid(timeArray, R);
		CaArray  = interp2(TT, RR, CaGrid, tt, rr, 'spline');
		CaSpline = spline(timeArray, CaArray);
		Ca       = @(t) ppval( CaSpline, t);

		DF = @(t, x) [  km1 * x(2) -  kp1 * Ca(t)        .* x(1); ... 
						km2 * x(3) - (kp2 * Ca(t) + km1) .* x(2) + kp1 * Ca(t) .* x(1); ...
						km3 * x(4) - (kp3 * Ca(t) + km2) .* x(3) + kp2 * Ca(t) .* x(2); ...
						km4 * x(5) - (kp4 * Ca(t) + km3) .* x(4) + kp3 * Ca(t) .* x(3); ...
						km5 * x(6) - (kp5 * Ca(t) + km4) .* x(5) + kp4 * Ca(t) .* x(4); ...
								   - (gamma       + km5) .* x(6) + kp5 * Ca(t) .* x(5) ];

		Jac = @(t, x) [              -  kp1 * Ca(t),         km1, 0, 0, 0, 0; ... 
                       kp1 * Ca(t),  - (kp2 * Ca(t) + km1),  km2, 0, 0, 0; ...
                    0, kp2 * Ca(t),  - (kp3 * Ca(t) + km2),  km3, 0, 0; ...
                 0, 0, kp3 * Ca(t),  - (kp4 * Ca(t) + km3),  km4, 0; ...
              0, 0, 0, kp4 * Ca(t),  - (kp5 * Ca(t) + km4),  km5; ...
           0, 0, 0, 0, kp5 * Ca(t),  - (gamma       + km5) ];

        [Tout, Xout] = odeSolve( DF, [0 totalTime], EQ, odeset('RelTol', ACC, 'Jacobian', Jac));
		Exo          = 1 - sum(Xout, 2); 
			
		ModelOutput = ModelOutput + fraction(Iter) * spline(Tout, Exo, DT);
    end  
end
