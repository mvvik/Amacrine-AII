function [Cost, Cost0] = CostFunction(GetModel, GetEquil, X, timeArray1, timeArray2, timeArray3, CaGrid1, CaGrid2, CaGrid3, DataIn, ErrorIn, rArray, DT)
%
%   Compute total model error vs model parameters: Cost = -log(Likelihood)
%   -------------------------------------------------------------------------
%   GetModel             = model exocytosis vs pulse duration function handle
%   GetEquil             = equilibrium of model Ca2+ sensor function handle
%   X                    = model parameters 
%   indKD_CP             = position of (KD, beta) parameters within X
%   [timeArray1 CaGRid1] = [Ca2+] vs time in 2mM  EGTA
%   [timeArray2 CaGRid2] = [Ca2+] vs time in 10mM EGTA
%   [timeArray3 CaGRid3] = [Ca2+] vs time in 1mM  BAPTA
%   DataIn               = capacitance vs pulse duration data, all conditions
%   ErrorIn              = capacitance vs pulse duration data, standard error
%   rArray               = distance array
%   DT                   = pulse duration array

    X        = SetParamBounds( X, 0 );  % Enforce parameter bounds
    EQ       = GetEquil(X);             % Get sensor equilibrium state
    Bnd      = EQ(end);                 % Equilibrium fully bound fraction will be penalized
    odeSolve = @ode15s;                 % Use stiff ODE solver
	ACC      = 1e-5;                    % Tolerance for the ODE solver

    Out1 = GetModel( EQ, CaGrid1, timeArray1, rArray, X, DT, odeSolve, ACC );  % Model exocytosis vs pulse duration for 2mM  EGTA
    Out2 = GetModel( EQ, CaGrid2, timeArray2, rArray, X, DT, odeSolve, ACC );  % Model exocytosis vs pulse duration for 10mM EGTA
    Out3 = GetModel( EQ, CaGrid3, timeArray3, rArray, X, DT, odeSolve, ACC );  % Model exocytosis vs pulse duration for 1mM  BAPTA
    
    DataOut = [Out1, Out2, Out3];                                                     % Combine model output across all buffering conditions
    scale   = sum(DataIn .* DataOut .* ErrorIn) / sum(DataOut.* DataOut .* ErrorIn);  % Scale model output to match experimental data
    Cost0   = sum( (DataIn - scale * DataOut).^2 .* ErrorIn);                         % Total model error = -2*log(Likelihood) +/- const
    Cost    = Cost0 * (1 + 20 * Bnd^2);                                               % Penalize fully bound sensor at background Ca2+ level
end
