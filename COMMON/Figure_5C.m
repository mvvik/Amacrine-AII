% -------------------------------------------------------------------------
%  Figure 5C: Capacitance vs pulse duration: experiment and model fits
% -------------------------------------------------------------------------

clear
if contains(pwd,'COMMON')
    fprintf('Do not execute scripts in this folder (switch to a Figure folder)\n');
    return;
end

addpath('../COMMON/');
ComputeSetup;
ResultsOut = CollectParForData(nPars, DataPrefix);

% -------------------------------------------------------------------------
% Plotting colors and styles
% -------------------------------------------------------------------------

clrBlue     = [0.6  0.6   1];
clrBlueLT   = [0.94 0.94  1];
clrOrange   = [1    0.6   0];
clrOrangeLT = [1    0.94  0.85];
clrE2       = [0    0     0];        % EGTA 2mM model curve
clrE10      = [0.7  0.3   1];        % EGTA 10mM model curve
clrB1       = [1    0.3   0.7];      % BAPTA 1mM model curve

ls = '--';  % Line style for model curves

% -------------------------------------------------------------------------
% Experimental data and statistics
% -------------------------------------------------------------------------

DTexp   = [0.5  1  5  10  20];   % Pulse durations

DataExp = [ ...
     6.328  24.537  49.24  53.739  62.329, ...   % EGTA 2mM
     9.706  28.063  54.52  57.003  69.097, ...   % EGTA 10mM
     0.043   6.951  30.436 37.228  37.281 ];     % BAPTA 1mM

SEM = [ ...
     3.672 12.000 25.645 24.668 28.762, ...
     6.834 12.064 19.755 25.471 53.280, ...
     2.509  3.840 21.645 27.532 14.239 ];

NN  = [10 10 16 15 37, 11 11 11 11 11, 10 7 9 8 8];  % Sample sizes
SEM = SEM ./ sqrt(NN);                               % Convert SEM to SE

% -------------------------------------------------------------------------
% Averaging EGTA data and preparing BAPTA data
% -------------------------------------------------------------------------

DataIn1 = 0.5 * (DataExp(1:5) + DataExp(6:10));        % Mean of EGTA 2 & 10 mM
SEM1    = 0.5 * sqrt(SEM(1:5).^2 + SEM(6:10).^2);      % Pooled SE
DataIn2 = DataExp(11:15);                              % BAPTA data
SEM2    = SEM(11:15);                                  % BAPTA SE

titles = {'EGTA (2mM)', 'EGTA (10mM)', 'BAPTA (1mM)'};

% -------------------------------------------------------------------------
% Compute model outputs
% -------------------------------------------------------------------------

fs     = 14;   % Axes   font size
lfs    = 15;   % Label  font size
tfs    = 18;   % Title  font size
lgndFs = 10;   % Legend font size

Vmin   = ResultsOut(1, end);     % Minimum cost
X      = ResultsOut(1, :);       % Best-fit parameter set

X      = SetParamBounds(X(1:nPars), 0);
EQ     = GetEquil(X);            % Steady-state sensor state
Bnd    = EQ(end);                % Final bound fraction

odeSolve = @ode15s;
ACC      = 1e-5;                 % Integration accuracy

Out1 = GetModel(EQ, CaGrid1, tArray1, rArray, X, DT, odeSolve, ACC);
Out2 = GetModel(EQ, CaGrid2, tArray2, rArray, X, DT, odeSolve, ACC);
Out3 = GetModel(EQ, CaGrid3, tArray3, rArray, X, DT, odeSolve, ACC);

% -------------------------------------------------------------------------
% Cost evaluation: weighted least squares
% -------------------------------------------------------------------------

DataOut = [Out1, Out2, Out3];

Range1 = 1 : nPoints;
Range2 = Range1 + 2*nPoints;

In1  = DataIn(Range1);
Err1 = 1 ./ sqrt(ErrorIn(Range1));

In3  = DataIn(Range2);
Err3 = 1 ./ sqrt(ErrorIn(Range2));

scale = sum(DataIn .* DataOut .* ErrorIn) / sum(DataOut.^2 .* ErrorIn);
cost  = sum((DataIn - scale * DataOut).^2 .* ErrorIn) * (1 + 20 * Bnd^2);

fprintf('Minimal model error (across all buffering conditions): %g \n', cost);

% -------------------------------------------------------------------------
% Plot: Experimental vs Model Data
% -------------------------------------------------------------------------

figure;

% EGTA: Fill region and plot error bars
fill([DT, fliplr(DT)], [In1 - Err1, fliplr(In1 + Err1)], clrBlueLT, 'EdgeColor', clrBlueLT); hold on;
errorbar(DTexp, DataIn1, SEM1, 'o', 'MarkerSize', 9, 'LineWidth', 2, 'Color', clrBlue);
plot(DT, In1, '-', 'LineWidth', 3, 'Color', clrBlue);
plot(DT, scale * Out1, 'LineStyle', ls, 'Color', clrE2, 'LineWidth', 2);
plot(DT, scale * Out2, 'LineStyle', ls, 'Color', clrE10, 'LineWidth', 2);

% BAPTA: Fill region and plot error bars
fill([DT, fliplr(DT)], [In3 - Err3, fliplr(In3 + Err3)], clrOrangeLT, 'EdgeColor', clrOrangeLT);
errorbar(DTexp, DataIn2, SEM2, 'o', 'MarkerSize', 9, 'LineWidth', 2, 'Color', clrOrange);
plot(DT, In3, '-', 'LineWidth', 3, 'Color', clrOrange);
plot(DT, scale * Out3, 'LineStyle', ls, 'LineWidth', 2, 'Color', clrB1);

% -------------------------------------------------------------------------
% Final Plot Formatting
% -------------------------------------------------------------------------

set(gca, 'FontSize', fs);
xlabel('Pulse duration (ms)', 'FontSize', lfs);
ylabel('\DeltaC_m (fF)', 'FontSize', lfs);
set(gca, 'XScale', 'log');
set(gca, 'Box', 'off');
set(gca, 'XTick', [0.5 1 2 5 10 20]);
axis tight;

legend( ...
    '', '', ...
    'EGTA 2mM & 10mM Data', ...
    'EGTA 2mM Model', ...
    'EGTA 10mM Model', ...
    '', '', ...
    'BAPTA 1mM Data', ...
    'BAPTA 1mM Model', ...
    'FontSize', lgndFs, 'Box', 'off', 'Location', 'NW' );

title('Best model fit:', 'FontSize', tfs);
fprintf(['Pars = [ ', fStr, '] \n'], X);


