%
% Global definitions: exocytosis model type, etc.
% -------------------------------------------------
%  1-pool cooperative exocytosis sensor parameters:
%     1     2       3     4    5 
%    R1  SCaKD1  SCakp1 coop gamma

addpath('../COMMON/');
ComputeCalciumPDE;

nCaSites   = 5;       % Number of Ca binding sites in exocytosis sensor (2..5)
nPars      = 5;       % Number of exocytosis model parameters
DataPrefix = ['Coop_1pool_CaSites', num2str(nCaSites), '_']; 

% --- Generate fprintf format string for parameter output:

fStr = '';
for ind = 1 : nPars
   fStr = [fStr, '%g '];
end

% --- Select the model integrator based on the number of sensor's Ca2+ binding sites:

GetModel = str2func( strcat('@GetModelOutput',         num2str(nCaSites)) );
GetEquil = str2func( strcat('@(X) getEquilibriumCoop', num2str(nCaSites), '(X(2), X(4))') );
CostFunc = @(X) CostFunction(GetModel, GetEquil, X, tArray1, tArray2, tArray3, CaGrid1, CaGrid2, CaGrid3, DataIn, ErrorIn, rArray, DT);

% -- Parameter windows to scan and initial error target before optimization:

switch nCaSites
    case 2
        Xmin = [ 0.035  0.5    0.4   0.003   1];
        Xmax = [ 0.050    4   0.99   0.50   10];
        maxCost = 135;
    case 3
        Xmin = [ 0.025   1     0.4   0.01   1];
        Xmax = [ 0.040   7     0.99   0.3  10];
        maxCost = 140;
    case 4
        Xmin = [ 0.020   3     0.4   0.02   1];
        Xmax = [ 0.030  11     0.99  0.50  15];
        maxCost = 150;
    case 5
        Xmin = [ 0.015   4     0.4   0.10   1];
        Xmax = [ 0.026  14     0.99  0.55  30];
        maxCost = 165;
end  

Xmin = real( SetParamBounds(Xmin, 1) );  % Invert bound-clipping to ensure the whole range is covered
Xmax = real( SetParamBounds(Xmax, 1) );  % Invert bound-clipping to ensure the whole range is covered


