function [ModelOutput, N2outDT, C2outDT] = Get_CaM_binding( CaGrid, timeArray, rArray, X, DT, odeSolve, ACC )

    totalTime = timeArray(end);

    [TT, RR] = meshgrid(timeArray, rArray);
    [tt, rr] = meshgrid(timeArray, X);
    CaArray  = interp2(TT, RR, CaGrid, tt, rr, 'spline');
    CaSpline = spline(timeArray, CaArray);
    Ca       = @(t) ppval( CaSpline, t);

    kpN1 = 0.1;   kpN2 = 0.15;  KN1 = 26.6;  KN2 = 6.6;
    kpC1 = 0.004; kpC2 = 0.01;  KC1 = 10;    KC2 = 0.93;
    
    kmN1 = kpN1 * KN1;
    kmN2 = kpN2 * KN2;
    kmC1 = kpC1 * KC1;
    kmC2 = kpC2 * KC2;

    EQN = Get_Equilibrium(KN1, KN2);
    EQC = Get_Equilibrium(KC1, KC2);

    DFN = @(t, x)  [    kmN1 * x(2) -  kpN1 * Ca(t)         .* x(1); ... 
					    kmN2 * x(3) - (kpN2 * Ca(t) + kmN1) .* x(2) + kpN1 * Ca(t) .* x(1); ...
                                  -                   kmN2  .* x(3) + kpN2 * Ca(t) .* x(2)     ];

    DFC = @(t, x)  [    kmC1 * x(2) -  kpC1 * Ca(t)         .* x(1); ... 
					    kmC2 * x(3) - (kpC2 * Ca(t) + kmC1) .* x(2) + kpC1 * Ca(t) .* x(1); ...
                                 -                    kmC2  .* x(3) + kpC2 * Ca(t) .* x(2)     ];

    JacN = @(t, x) [                    -  kpN1 * Ca(t),          kmN1, 0; ... 
				         kpN1 * Ca(t),  - (kpN2 * Ca(t) + kmN1),  kmN2;    ...
	                  0, kpN2 * Ca(t),  -               + kmN2                 ];


	JacC = @(t, x) [                    -  kpC1 * Ca(t),          kmC1, 0; ... 
				         kpC1 * Ca(t),  - (kpC2 * Ca(t) + kmC1),  kmC2;    ...
	                  0, kpC2 * Ca(t),  -               + kmC2                 ];


    [TNout, Nout] = odeSolve( DFN, [0 totalTime], EQN, odeset('RelTol', ACC, 'Jacobian', JacN));
    [TCout, Cout] = odeSolve( DFC, [0 totalTime], EQC, odeset('RelTol', ACC, 'Jacobian', JacC));

    N2outDT = spline(TNout, Nout(:,end), DT);
    C2outDT = spline(TCout, Cout(:,end), DT);
    ModelOutput = N2outDT .* C2outDT;
end

