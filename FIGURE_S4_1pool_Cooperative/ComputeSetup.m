%
% Global definitions: exocytosis model type, etc.
% -------------------------------------------------
%  1-pool cooperative exocytosis sensor parameters:
%     1     2       3     4    5 
%    R1  SCaKD1  SCakp1 coop gamma

addpath('../COMMON/');
ComputeCalciumPDE;

nPars      = 5;       % Number of exocytosis model parameters
nCaSites   = 5;       % Number of Ca binding sites in exocytosis sensor (2..5)
DataPrefix = ['Coop_1pool_CaSites', num2str(nCaSites), '_']; 

% --- Generate exocytosis data curve that the model will attempt fitting:

DT                = logspace(log10(0.5), log10(totalTime), nPoints); 
[DataIn, ErrorIn] = GenerateExocytosisData( DT, 0);

% --- Read computed [Ca2+] vs distance and time, for each buffering condition:

[  ~,    CaGrid1, tArray1] = ReadDataVsTime(fileEGTA2mM );
[  ~,    CaGrid2, tArray2] = ReadDataVsTime(fileEGTA10mM);
[rArray, CaGrid3, tArray3] = ReadDataVsTime(fileBAPTA1mM);

% --- Generate fprintf format string for parameter output:

fStr = '';
for ind = 1 : nPars
   fStr = [fStr, '%g '];
end

% --- Select the model integrator based on the number of sensor's Ca2+ binding sites:

GetModel = str2func( strcat('@GetModelOutput',         num2str(nCaSites)) );
GetEquil = str2func( strcat('@(X) getEquilibriumCoop', num2str(nCaSites), '(X(2), X(4))') );
CostFunc = @(X) CostFunction(GetModel, GetEquil, X, tArray1, tArray2, tArray3, CaGrid1, CaGrid2, CaGrid3, DataIn, ErrorIn, rArray, DT);

% --- Parameter windows to scan:

switch nCaSites
    case 2
        Xmin = [ 0.035  0.5    0.4   0.003   1];
        Xmax = [ 0.050    4   0.99   0.50   10];
    case 3
        Xmin = [ 0.025   1     0.4   0.01   1];
        Xmax = [ 0.040   7     0.99   0.3  10];
    case 4
        Xmin = [ 0.020   3     0.4   0.02   1];
        Xmax = [ 0.030  11     0.99  0.50  15];
    case 5
        Xmin = [ 0.015   4     0.4   0.10   1];
        Xmax = [ 0.026  14     0.99  0.55  30];
end  

Xmin = real( SetParamBounds(Xmin, 1) );  % Invert bound-clipping to ensure the whole range is covered
Xmax = real( SetParamBounds(Xmax, 1) );  % Invert bound-clipping to ensure the whole range is covered

maxCost = 180;

