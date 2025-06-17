clear;		
%ComputeCalciumPDE;
nCaSites = 5;
DefineGlobalParameters;

nTrials1 = 100;
nTrials2 = 150;

OutName       = GetNextFileName( prefixN ); 
[~, HostName] = system('hostname');
save(OutName);

[  ~,    CaGrid1, tArray1] = ReadDataVsTime(nameEGTA2mM );
[  ~,    CaGrid2, tArray2] = ReadDataVsTime(nameEGTA10mM);
[rArray, CaGrid3, tArray3] = ReadDataVsTime(nameBAPTA1mM);

%   Parameters: distances, sensor binding kinetics, and pool fraction:
%      1      2      3        4       5    6        7   
%   Dist-1 Dist-2 SensorKD SensorKp coop gamma pool-fraction

switch nCaSites
    case 5
        Xmin = [ 0.020  0.025   0.5  0.20  0.05  1.0  0.50];
        Xmax = [ 0.035  0.060  11.0  0.99  0.60  50   0.70];
        GetModel = @GetModelOutput5;
        GetEquil = @getEquilibriumCoop5;
    case 4
        Xmin = [ 0.023  0.035  0.3  0.45  0.04  2.0  0.48];
        Xmax = [ 0.035  0.065  4.0  0.99  0.60  12   0.68];
        GetModel = @GetModelOutput4;
        GetEquil = @getEquilibriumCoop4;
    case 3
        Xmin = [ 0.028  0.045  0.2   0.45  0.01  2.0  0.48];
        Xmax = [ 0.040  0.080  2.5   0.99  0.50  7.0  0.70];
        GetModel = @GetModelOutput3;
        GetEquil = @getEquilibriumCoop3;
    case 2
        Xmin = [ 0.035  0.065  0.3   0.55  0.01  2.0  0.50];
        Xmax = [ 0.050  0.010  1.4   0.99  0.40  4.0  0.75];
        GetModel = @GetModelOutput2;
        GetEquil = @getEquilibriumCoop2;
end

maxCost    = 150;
nTrials    = nTrials1 * nTrials2;
ResultsOut = zeros(nTrials1, nTrials2, nPars+1);

Xmin = real( SetParamBounds(Xmin, Rpatch, 1) );
Xmax = real( SetParamBounds(Xmax, Rpatch, 1) );  Xmax(4) = 2.0;

options = optimset('Display', 'off', 'TolFun', 3e-3, 'TolX', 1e-4, 'MaxFunEvals', 1500);
%rs      = RandStream.create('mrg32k3a', 'NumStreams', nTrials1, 'Seed', 'shuffle', 'CellOutput', true);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parfor n1 = 1 : nTrials1
   
    %str = RandStream('mt19937ar','Seed', 0   )
    Cost0 = 2*maxCost;
    str   = RandStream('mrg32k3a', 'Seed', 3*n1 + 1);
    RandStream.setGlobalStream(str);

    while Cost0 > maxCost
        P0    = exp( unifrnd( log(Xmin), log(Xmax) ) );
        Cost0 = CostFunction(GetModel, GetEquil, P0, tArray1, tArray2, tArray3, CaGrid1, CaGrid2, CaGrid3, DataIn, ErrorIn, rArray, DT);
    end
    fprintf("Accepted [%g %g %g %g %g %g %g] Cost0 = %g \n", P0, Cost0);

    Results = zeros(nTrials2, nPars+1);
    n2      = 0;
    fff     = 0.025;

    while n2 < nTrials2

        %fff = 0.005 + 0.05 * (1 - n2 / nTrials2);
        
        for k = 1 : 100
            Pars = unifrnd(P0*(1-fff), P0*(1+fff) );
            Cost = CostFunction(GetModel, GetEquil, Pars, tArray1, tArray2, tArray3, CaGrid1, CaGrid2, CaGrid3, DataIn, ErrorIn, rArray, DT);
            test = exp((Cost0 - Cost)/2);
            if rand < test
               break;
            end
        end
        
        if Cost < maxCost
            n2 = n2 + 1;
            Results(n2, :) = [Pars, Cost];
        end
        P0 = Pars; Cost0 = Cost;
    end
  
    [C, I] = min(Results(:,end));
    fprintf("Min: [%.3g %.3g %.3g %.3g %.3g %.3g %.3g]: Cost = %g \n", Results(I,:));
    ResultsOut(n1, :, :) = Results;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Results = reshape(ResultsOut, nTrials, nPars + 1);
inds = find(ResultsOut(:,end) < maxCost);
ResultsOut = ResultsOut(inds, :);
SaveParForData;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C0      = min(ResultsOut(:, end));
C1      = C0 + 8.5;
inds    = find(ResultsOut(:,end) < C1);
ResultsOut = ResultsOut(inds, :);

[~, ind] = min(ResultsOut(:,end));
xBest    = ResultsOut(ind, 1:nPars);

nBins    = 50; 
nTotal   = nPars * nBins;
Results  = zeros(nPars, nBins, nPars + 1);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for indPar = 1 : 7

    parfor indBin = 1 : nBins

        valPar0 = Xmin(indPar) + (indBin - 1) / nBins * (Xmax(indPar) - Xmin(indPar));
        valPar1 = Xmin(indPar) +  indBin      / nBins * (Xmax(indPar) - Xmin(indPar));
        valPar  = 0.5 * (valPar0 + valPar1);

        if indPar == 6
            valPar0 = exp(valPar0);
            valPar1 = exp(valPar1);
            valPar  = exp( valPar);
        end

        fprintf("Par = %d, val = %g .. %g \n", indPar, valPar0, valPar1);
        
        inds1 = find(ResultsOut(:, indPar) >= valPar0);
        inds2 = [];

        if numel(inds1) > 0
            inds2 = find(ResultsOut(inds1, indPar) < valPar1);
        end
        if numel(inds2) > 0
            RES       = ResultsOut(inds1(inds2), :);
            [xx, ind] = min(RES(:,end));
            X0        = RES(ind, 1:nPars);
        else
            [C, ind]  = min( (ResultsOut(:, indPar) - valPar).^2 );
            X0        = ResultsOut(ind, 1:nPars);
        end

        cost0 = CostFunction(GetModel, GetEquil, X0, tArray1, tArray2, tArray3, CaGrid1, CaGrid2, CaGrid3, DataIn, ErrorIn, rArray, DT);

        fprintf(['> Starting cost0 = %g;  Pars: ', fStr, '\n'], cost0, X0);
        X  = removeValue(X0,    indPar);
        [X,  cost ] = fminsearch( @(X) CostFunction(GetModel, GetEquil, insertValue(X, indPar, valPar), tArray1, tArray2, tArray3, CaGrid1, CaGrid2, CaGrid3, DataIn, ErrorIn, rArray, DT), X,  options );
        Y = insertValue(X, indPar, valPar);

        Results(indPar, indBin, :) = [Y, cost];
        fprintf(['  finished %g --> %g \n  Pars: ', fStr, '\n\n'], cost0, cost, Y);
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ResultsOld = Results;
Results    = reshape( Results, nTotal, nPars + 1);
SaveParForData;


