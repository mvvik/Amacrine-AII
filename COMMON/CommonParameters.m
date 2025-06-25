%
% Ca2+ simulation and other common global parameters, except for exocytosis sensor parameters
%
CaBgr      = 0.05;                  % Background Ca concentration, uM
Rpatch     = 0.25;                  % Simulation domain radius, um
totalTime  = 20;                    % Total simulation time, ms
Kappa      = 50;                    % Endogenous buffering capacity
BKD        = 1;                     % Endogenous buffer Ca2+ affinity, uM
BDiff      = 0.05;                  % Buffer mobility (um^2/ms)
Bkplus     = 0.4;                   % Endogenous buffer Ca2+ binding rate, 1/(uM ms)
ICa        = 0.2;                   % Ca2+ current (times 2 due to reflection symmetry), pA
CaTau0     = 10;                    % Total free Ca2+ clearance rate, 1/ms
CaTau      = CaTau0 / (1 + Kappa);  % Effective extrusion rate with buffer
nGrid      = 100;                   % Number of spatial nodes
nPoints    = 50;                    % Number of time points in simulation to fit to data
LLR_CutOff = 8.0;                   % Log likelihood ratio cut-off
