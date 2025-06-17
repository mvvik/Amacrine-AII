
%   Parameters/X: distances and sensor parameters, plus pool fractions:
%      1      2       3         4        5     6        7   
%   Dist-1 Dist-2  SensorKD  SensorKp  coop  gamma  pool-fraction

clear;
nCaSites = 5;
GetModel = @GetModelOutputShorter5;
GetEquil = @getEquilibriumCoop5;
CollectParForData;

[     ~, CaGrid1, tArray1] = ReadDataVsTime(nameEGTA2mM);
[     ~, CaGrid2, tArray2] = ReadDataVsTime(nameEGTA10mM);
[rArray, CaGrid3, tArray3] = ReadDataVsTime(nameBAPTA1mM);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clrBlue     = [0.6 0.6     1];
clrBlueLT   = [0.94 0.94   1];
clrOrange   = [1   0.6     0];
clrOrangeLT = [1   0.94  0.85];
clrE2       = [0    0      0];
clrE10      = [0.7  0.3    1];
clrB1       = [1    0.3  0.7];

ls = '--';

DTexp   = [ 0.5      1     5       10     20];

DataExp  = [6.328  24.537  49.24  53.739  62.329...
           9.706  28.063  54.52  57.003  69.097...
           0.043   6.951  30.436 37.228	 37.281];

SEM     = [3.672  12.000  25.645 24.668  28.762...
           6.834  12.064  19.755 25.471	 53.280...
           2.509   3.840  21.645 27.532  14.239];
       
NN      = [10 10 16 15 37 ...
           11 11 11 11 11 ...
           10  7  9  8 8 ]; 

SEM     = SEM ./ sqrt(NN);       
       
DataIn1 = 0.5 * ( DataExp(1 : 5) + DataExp(6 : 10) );
SEM1    = 0.5 * sqrt( SEM(1 : 5).^2 + SEM(6 : 10).^2); 
DataIn2 = DataExp( 11 : 15 );
SEM2    =     SEM( 11 : 15 );

titles = { 'EGTA (2mM)', 'EGTA (10mM)', 'BAPTA (1mM)'};

fs     = 14;
lfs    = 15;
lgndFs = 10;

figure; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for kk = 2

    Vmin = ResultsOut(1, end);

    switch kk
        case 1
            X = ResultsOut(1, :);
        case 2
            ind = find(ResultsOut(:,end) > Vmin + 2, 1);
            X   = ResultsOut(ind, :);
    end

    X   = SetParamBounds( X(1:nPars), Rpatch, 0 );
    EQ  = GetEquil( X(3), X(5) );
    Bnd = EQ(end);   

    odeSolve = @ode15s;
	ACC      = 1e-5;

    Out1 = GetModel( EQ, CaGrid1, tArray1, rArray, X, DT, odeSolve, ACC );
    Out2 = GetModel( EQ, CaGrid2, tArray2, rArray, X, DT, odeSolve, ACC );
    Out3 = GetModel( EQ, CaGrid3, tArray3, rArray, X, DT, odeSolve, ACC );
    
    DataOut = [Out1, Out2, Out3];
    Range1 = 1 : nPoints;
    Range2 = Range1 + 2*nPoints;
    In1 = DataIn(Range1);   Err1 = 1./sqrt(ErrorIn(Range1));
    In3 = DataIn(Range2);   Err3 = 1./sqrt(ErrorIn(Range2));
   
    scale = sum(DataIn .* DataOut .* ErrorIn) / sum(DataOut.* DataOut .* ErrorIn);
    cost  = sum( (DataIn - scale * DataOut).^2 .* ErrorIn) * (1 + 20 * Bnd^2);
    
    fprintf(' Cost: %g \n', cost ); 
 
    fill([DT, fliplr(DT)], [In1-Err1, fliplr(In1+Err1)],  clrBlueLT, 'EdgeColor', clrBlueLT); hold on;
    errorbar(DTexp, DataIn1, SEM1, 'o', 'MarkerSize', 9, 'linewidth', 2, 'color', clrBlue); 
    plot(DT,       In1,  '-', 'linewidth', 3,  'color', clrBlue); hold on;
    plot(DT, scale*Out1, 'linestyle', ls, 'color', clrE2,  'linewidth', 2); 
    plot(DT, scale*Out2, 'linestyle', ls, 'color', clrE10, 'linewidth', 2);
    
    fill([DT, fliplr(DT)], [In3-Err3, fliplr(In3+Err3)],  clrOrangeLT, 'EdgeColor', clrOrangeLT); hold on;
    errorbar(DTexp, DataIn2, SEM2, 'o', 'MarkerSize', 9, 'linewidth', 2, 'color', clrOrange); 
    plot(DT,       In3,  '-', 'linewidth', 3, 'color', clrOrange); 
    plot(DT, scale*Out3, 'linestyle', ls, 'linewidth', 2, 'color', clrB1); 
    
    set( gca,                     'fontsize', fs );
    xlabel('Pulse duration (ms)', 'fontsize', lfs);
    ylabel('\DeltaC_m (fF)',      'fontsize', lfs);
    legend('', '', 'EGTA 2mM & 10mM Data', 'EGTA 2mM Model', 'EGTA 10mM Model', '', '', 'BAPTA 1mM Data', 'BAPTA 1mM Model', 'fontsize', lgndFs, 'box', 'off');
    set(gca, 'xscale', 'log');
    set(gca, 'box', 'off');
    set(gca,'xtick', [1 2 5 10]);
    axis tight;
    fprintf(['Pars = [ ', fStr, '] \n'], X );
   
    set(gca, 'xtick', [0.5 1 2 5 10 20]);
    title(sprintf('Cost = %.3g, Bnd = %.3g', cost, Bnd));
end
