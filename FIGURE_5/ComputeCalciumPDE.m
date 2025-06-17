DefineGlobalParameters;

tSteps    = totalTime * 100;

Pars = [ Rpatch, Kappa, BKD, Bkplus, CaTau, ICa, nGrid, totalTime, tSteps ];

[X, Y] = system( [ prog, ' ', num2str([2e3, 0, Pars]), ' ', nameEGTA2mM  ]);  % EGTA   2mM
[X, Y] = system( [ prog, ' ', num2str([1e4, 0, Pars]), ' ', nameEGTA10mM ]);  % EGTA   10mM
[X, Y] = system( [ prog, ' ', num2str([0, 1e3, Pars]), ' ', nameBAPTA1mM ]);  % BAPTA  1mM
