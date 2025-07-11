%
% Global definitions: exocytosis model type, etc.
% -------------------------------------------------
% 1-pool non-cooperative exocytosis sensor parameters:
%     1     2       3     4    
%    R1  SCaKD1  SCakp1 gamma   

addpath('../COMMON/');
ComputeCalciumPDE;

nCaSites   = 5;       % Number of Ca binding sites in exocytosis sensor (2..5)
nPars      = 4;       % Number of exocytosis model parameters
DataPrefix = ['Noncoop_1pool_CaSites', num2str(nCaSites), '_']; 

% --- Generate fprintf format string for parameter output:

fStr = '';
for ind = 1 : nPars
   fStr = [fStr, '%g '];
end

% --- Select the model integrator based on the number of sensor's Ca2+ binding sites:

GetModel = str2func( strcat('@GetModelOutput',         num2str(nCaSites)) );
GetEquil = str2func( strcat('@(X) getEquilibriumCoop', num2str(nCaSites), '(X(2), 1)') );
CostFunc = @(X) CostFunction(GetModel, GetEquil, X, tArray1, tArray2, tArray3, CaGrid1, CaGrid2, CaGrid3, DataIn, ErrorIn, rArray, DT);

% -- Parameter windows to scan and initial error target before optimization:

switch nCaSites
        case 2
            Xmin = [  0.04,  0.4,  0.45,    1];
	    Xmax = [  0.07,    4,  0.99,   10];
            maxCost = 180;
        case 3
            Xmin = [  0.030,  1,   0.4,     1];
	    Xmax = [  0.055,  6,  0.99,   200];
            maxCost = 165;
        case 4
            Xmin = [  0.025,  1,   0.3,    2.5];
	    Xmax = [  0.045,  7,  0.99,   9000];
            maxCost = 155;
        case 5
            Xmin = [  0.023,   1,   0.3,  3.5];
	    Xmax = [  0.033,   7,  0.99,  9000];
            maxCost = 150;
end

Xmin = real( SetParamBounds(Xmin, 1) );  % Invert bound-clipping to ensure the whole range is covered
Xmax = real( SetParamBounds(Xmax, 1) );  % Invert bound-clipping to ensure the whole range is covered



