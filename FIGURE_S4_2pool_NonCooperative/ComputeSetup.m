%
% Global definitions: exocytosis model type, etc.
% -------------------------------------------------
% 2-pool non-cooperative exocytosis sensor parameters:
%      1      2      3        4       5    6       
%   Dist-1 Dist-2 SensorKD SensorKp gamma pool-fraction

addpath('../COMMON/');
ComputeCalciumPDE;

nCaSites   = 5;       % Number of Ca binding sites in exocytosis sensor (2..5)
nPars      = 6;       % Number of exocytosis model parameters
DataPrefix = ['Noncoop_2pool_CaSites', num2str(nCaSites), '_']; 

% --- Generate fprintf format string for parameter output:

fStr = '';
for ind = 1 : nPars
   fStr = [fStr, '%g '];
end

% --- Select the model integrator based on the number of sensor's Ca2+ binding sites:

GetModel = str2func( strcat('@GetModelOutput',         num2str(nCaSites)) );
GetEquil = str2func( strcat('@(X) getEquilibriumCoop', num2str(nCaSites), '(X(3), 1)') );
CostFunc = @(X) CostFunction(GetModel, GetEquil, X, tArray1, tArray2, tArray3, CaGrid1, CaGrid2, CaGrid3, DataIn, ErrorIn, rArray, DT);

% -- Parameter windows to scan and initial error target before optimization:

switch nCaSites
        case 5
            Xmin = [ 0.026  0.046  0.5   0.5    1   0.54];
            Xmax = [ 0.037  0.066  2.0  0.99   30   0.65];
            maxCost = 70;
        case 4
            Xmin = [ 0.031  0.050  0.6   0.4    1   0.48];
            Xmax = [ 0.040  0.070  1.4  0.99   20   0.70];
            maxCost = 75;
        case 3
            Xmin = [ 0.035  0.060  0.6  0.50   1  0.48];
            Xmax = [ 0.045  0.082  1.5  0.99  20  0.70];
            maxCost = 85;
        case 2
            Xmin = [ 0.028  0.045  0.55  0.55    1   0.48];
            Xmax = [ 0.040  0.065  1.60  0.99   20   0.70];
            maxCost = 110;
end

Xmin = real( SetParamBounds(Xmin, 1) );  % Invert bound-clipping to ensure the whole range is covered
Xmax = real( SetParamBounds(Xmax, 1) );  % Invert bound-clipping to ensure the whole range is covered

