%
% Global definitions: exocytosis model type, etc.
% -----------------------------------------------------------
% 2-pool cooperative exocytosis sensor parameters:
%      1      2      3        4       5    6        7   
%   Dist-1 Dist-2 SensorKD SensorKp coop gamma pool-fraction

addpath('../COMMON/');
ComputeCalciumPDE;    % Compute [Ca2+] vs time and space

nCaSites   = 5;       % Number of Ca binding sites in exocytosis sensor (2..5)
nPars      = 7;       % Number of exocytosis model parameters
DataPrefix = ['Coop_2pool_CaSites', num2str(nCaSites), '_']; 

% --- Generate fprintf format string for parameter output:

fStr = '';
for ind = 1 : nPars
   fStr = [fStr, '%g '];
end

% --- Select the model integrator based on the number of sensor's Ca2+ binding sites:

GetModel = str2func( strcat('@GetModelOutput',         num2str(nCaSites)) );
GetEquil = str2func( strcat('@(X) getEquilibriumCoop', num2str(nCaSites), '(X(3), X(5))') );
CostFunc = @(X) CostFunction(GetModel, GetEquil, X, tArray1, tArray2, tArray3, CaGrid1, CaGrid2, CaGrid3, DataIn, ErrorIn, rArray, DT);

% -- Parameter windows to scan and initial model error before optimization:

switch nCaSites
    case 5
        Xmin = [ 0.018  0.025  0.1  0.30  0.05  1.0  0.45];
        Xmax = [ 0.035  0.060  8.0  0.99  0.90  80   0.75];
        maxCost  = 60;
    case 4
        Xmin = [ 0.023  0.035  0.3  0.45  0.04  2.0  0.48];
        Xmax = [ 0.035  0.065  4.0  0.99  0.60  12   0.68];
        maxCost  = 65;
    case 3
        Xmin = [ 0.028  0.045  0.2   0.45  0.01  2.0  0.48];
        Xmax = [ 0.040  0.080  2.5   0.99  0.50  7.0  0.70];
        maxCost  = 70;
    case 2
        Xmin = [ 0.035  0.065  0.3   0.55  0.01  2.0  0.50];
        Xmax = [ 0.048  0.092  1.3   0.99  0.30  4.0  0.75];
        maxCost  = 95;
end

Xmin = real( SetParamBounds(Xmin, 1) );  % Invert bound-clipping to ensure the whole range is covered
Xmax = real( SetParamBounds(Xmax, 1) );  % Invert bound-clipping to ensure the whole range is covered
