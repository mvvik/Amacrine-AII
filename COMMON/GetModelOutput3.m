function ModelOutput = GetModelOutput3( EQ, CaGrid, timeArray, rArray, X, DT, odeSolve, ACC )
%
%  Generate capacitance vs pulse duration data for a sensor with 3 binding sites
%  ------------------------------------------------------------------------------
%  EQ        = equilibrium of all sensor states at background Ca2+ level
%  CaGrid    = [Ca2+] vs time and distance
%  timeArray = temporal grid
%  rArray    = spatial grid
%  X         = model parameters (from 4 to 7, depending on the model)
%  DT        = Pulse duration array 
%  odeSolve  = ODE solver to use (ode15s is best)
%  ACC       = ODE solver tolerance (should be small for ode15s, e.g. 1e-5)
    
    totalTime   = DT(end);
	[TT, RR]    = meshgrid(timeArray, rArray);
			
    fraction  = 1;                    % Default value (single-pool model)
    nPars     = numel(X);             % No. of parameters --> model type and parameter encoding  
    nSites    = 1 + floor(nPars / 6); % Distinguish 1-pool model vs 2-pool model
    coopFlag  = bitand(nPars, 1);     % Distinguish cooperative vs. non-cooperative model
	
    if nSites == 2, fraction = [ X(nPars),  1 - X(nPars) ]; end  % Pool fractions in 2-pool model
    KD    = X(nSites + 1);                          % Sensor Ca2+ affinity
    kp    = X(nSites + 2);                          % Sensor Ca2+ binding rate
    coop  = X(nSites + 3) * coopFlag + ~coopFlag;   % Sensor Ca2+ binding cooperativity ("beta")
    gamma = X(nSites + 3 + coopFlag);               % Sensor final transition
    km    = KD * kp;                                % Sensor unbinding rate
    kp1   = 3*kp;     kp2 = 2*kp;  kp3 = kp         % Binding rates of individual stages
    km1   = km/coop;  km2 = 2*km;  km3 = 3*km*coop; % Unbinding rates of individual stages

	ModelOutput = zeros( size( DT ) );

    for Iter = 1 : nSites
	
		R        = X( Iter );
		[tt, rr] = meshgrid(timeArray, R);
		CaArray  = interp2(TT, RR, CaGrid, tt, rr, 'spline');
		CaSpline = spline(timeArray, CaArray);
		Ca       = @(t) ppval( CaSpline, t);
				
		DF = @(t, x) [  km1 * x(2) -  kp1 * Ca(t)        .* x(1); ... 
						km2 * x(3) - (kp2 * Ca(t) + km1) .* x(2) + kp1 * Ca(t) .* x(1); ...
						km3 * x(4) - (kp3 * Ca(t) + km2) .* x(3) + kp2 * Ca(t) .* x(2); ...
								   - (gamma       + km3) .* x(4) + kp3 * Ca(t) .* x(3) ];

		Jac = @(t, x) [          -  kp1 * Ca(t),         km1, 0, 0; ... 
				   kp1 * Ca(t),  - (kp2 * Ca(t) + km1),  km2, 0; ...
				0, kp2 * Ca(t),  - (kp3 * Ca(t) + km2),  km3; ...
	         0, 0, kp3 * Ca(t),  - (gamma       + km3) ];

        [Tout, Xout] = odeSolve( DF, [0 totalTime], EQ, odeset('RelTol', ACC, 'Jacobian', Jac));
        Exo          = 1 - sum(Xout, 2); 
		ModelOutput  = ModelOutput + fraction(Iter) * spline(Tout, Exo, DT);
    end  
end

