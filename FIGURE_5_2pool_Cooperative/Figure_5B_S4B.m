function Figure_5B_S4B()

% -------------------------------------------------------------------------
%   Figure 5B/S4B: Parameter profile plots (max likelihood profile) 
% -------------------------------------------------------------------------%

ComputeSetup;
ResultsOut = CollectParForData(nPars, DataPrefix);
if numel(ResultsOut) == 0
    fprintf("No data yet for model %s: run Likelihood_Profile first\n", DataPrefix);
    return;
end

%% Parameter Definitions
labels = { 'R_X', 'R_Y', 'K_D', 'k_{on}', '\beta', '\gamma', 'Pool_X', 'Pool_Y'};
Units  = { 'nm',  'nm',  '\muM', '(μM ms)^{-1}', '', 'ms^{-1}', '%', '%' }; 
Order  = [1 3 7 5, 2 4 8 6];    % Reordering of parameters for plotting

nBins  = 30;
clr    = [0.6, 0.0, 0.6];  % Fill color for confidence intervals
nRows  = 2;
nCols  = 4;
tfs    = 12;

indFrac  = nPars;         % Pool fraction parameter index (pre-augmentation)
indGamma = nPars + 1;     % Gamma parameter index (after augmentation)

%% Clamp parameter bounds and filter by error threshold
ResultsOut(:, 1:nPars) = SetParamBounds(ResultsOut(:, 1:nPars), 0);

[Vmin, Imin] = min(ResultsOut(:,end));
ParMin       = ResultsOut(Imin, 1:nPars);

fprintf('Minimal model exocytosis error (across all EGTA and BAPTA conditions): %g \n', Vmin);
fprintf([' Best-fit parameters = ', fStr, '\n'], ParMin);

MaxError           = Vmin + LLR_CutOff + 1;
ResultsOut         = ResultsOut(ResultsOut(:,end) <= MaxError, :);
ResultsOut(:, end) = ResultsOut(:,end) - Vmin;
nData              = size(ResultsOut, 1);
fprintf('Filtered set: N = %d\n', nData);

%% Rescale parameters for plotting
ResultsOut(:,1:2)     = ResultsOut(:,1:2)   * 1000;   % Convert distances to nm
ResultsOut(:,nPars)   = ResultsOut(:,nPars) * 100;    % Convert pool fraction to %

%% Add complementary pool fraction as a new parameter
nPars = nPars + 1;
ResultsNew              = zeros(nData, nPars + 1);
ResultsNew(:, 1:nPars ) = [ResultsOut(:, 1:nPars-1), 100 - ResultsOut(:, nPars-1)];
ResultsNew(:, end  )    = ResultsOut(:, end);
ResultsOut              = ResultsNew;
ParMin                  = ResultsOut(1, 1:nPars);

%% Swap R_X and R_Y values to enforce R_X < R_Y
for jj = 1:nData
    if ResultsOut(jj, 2) < ResultsOut(jj, 1)
        ResultsOut(jj, [1,2])       = ResultsOut(jj, [2,1]);
        ResultsOut(jj, [nPars-1,nPars]) = ResultsOut(jj, [nPars, nPars-1]);
    end
end

%% Reorder parameters for plotting
for jj = 1:nData
    ResultsOut(jj, 1:nPars) = ResultsOut(jj, Order); 
end

%% Plot confidence profiles for each parameter
figure;
for indPar = 1:nPars
    Xmin = min(ResultsOut(:, indPar));
    Xmax = max(ResultsOut(:, indPar));
    if indPar == indGamma
        Xmin = log(Xmin); Xmax = log(Xmax);
    end

    valPars = zeros(1, nBins);
    valErrs = zeros(1, nBins);

    for indBin = 1:nBins
        valPar0 = Xmin + (indBin - 1) / nBins * (Xmax - Xmin);
        valPar1 = Xmin + indBin       / nBins * (Xmax - Xmin);
        valPar  = 0.5 * (valPar0 + valPar1);

        if indPar == indGamma
            valPar0 = exp(valPar0); valPar1 = exp(valPar1); valPar = exp(valPar);
        end

        valPars(indBin) = valPar;
        inds = find(ResultsOut(:, indPar) >= valPar0 & ResultsOut(:, indPar) < valPar1);

        if ~isempty(inds)
            [~, idx] = min(ResultsOut(inds,end));
            valErrs(indBin) = ResultsOut(inds(idx), end);
        else
            [~, idx] = min((ResultsOut(:, indPar) - valPar).^2);
            valErrs(indBin) = ResultsOut(idx, end);
        end
    end

    if indPar == indGamma
        Xmin = exp(Xmin); Xmax = exp(Xmax);
    end

    subplot(nRows, nCols, indPar); hold on;
    fill([valPars(1), valPars, valPars(end)], [LLR_CutOff, valErrs, LLR_CutOff], clr, 'EdgeColor', clr, 'LineWidth', 1);
    axis([Xmin Xmax 0 LLR_CutOff]);
end

%% Overlay stratified samples by R_X (distance)
for mm = 0 : 2 : 26
   clrG = (0.9 - mm/40);  clrR  = 1 - clrG;  clr = [ clrR, clrG, 0 ];

   val0  = 19 + mm/2;
   val1  = 20 + mm/2;
   val2  = 21 + mm/2;
   inds1 = find(ResultsOut(:, 1) >= val0 & ResultsOut(:, 1) < val1 );
   inds2 = find(ResultsOut(:, 1) >= val1 & ResultsOut(:, 1) < val2 );
   [~, ind1] = min(ResultsOut(inds1,end));
   [~, ind2] = min(ResultsOut(inds2,end));
   Y1 = ResultsOut(inds1(ind1), end);
   Y2 = ResultsOut(inds2(ind2), end);
   X1 = ResultsOut(inds1(ind1), 1);
   X2 = ResultsOut(inds2(ind2), 1);
   Z1 = ResultsOut(inds1(ind1), 5);
   Z2 = ResultsOut(inds2(ind2), 5);

   subplot(nRows, nCols, 1); hold on;
   fill([X1 X1 X2 X2], [LLR_CutOff, Y1, Y2, LLR_CutOff], clr, 'EdgeColor', clr, 'LineWidth', 1);

   subplot(nRows, nCols, 5); hold on;
   fill([Z1 Z1 Z2 Z2], [LLR_CutOff, Y1, Y2, LLR_CutOff], clr, 'EdgeColor', clr, 'LineWidth', 1);
end

%% Overlay stratified samples by gamma
for mm = 1 : 18
   clrB = (0.9 - mm/20);  clrR  = 1 - clrB;  clr = [ clrR, 0, clrB ];

   val0  = 0.38 + 0.6*(mm - 1)/19;
   val1  = 0.38 + 0.6*(mm - 0)/19;
   val2  = 0.38 + 0.6*(mm + 1)/19;

   inds1 = find(ResultsOut(:, 6) >= val0 & ResultsOut(:, 6) < val1 );
   inds2 = find(ResultsOut(:, 6) >= val1 & ResultsOut(:, 6) < val2 );
   [~, ind1] = min(ResultsOut(inds1,end));
   [~, ind2] = min(ResultsOut(inds2,end));
   Y1 = ResultsOut(inds1(ind1), end);
   Y2 = ResultsOut(inds2(ind2), end);
   X1 = ResultsOut(inds1(ind1), 6);
   X2 = ResultsOut(inds2(ind2), 6);
   Z1 = ResultsOut(inds1(ind1), 8);
   Z2 = ResultsOut(inds2(ind2), 8);

   fprintf('mm = %d\n', mm);

   subplot(nRows, nCols, 6); hold on;
   fill([X1 X1 X2 X2], [LLR_CutOff, Y1, Y2, LLR_CutOff], clr, 'EdgeColor', clr, 'LineWidth', 1);

   subplot(nRows, nCols, 8); hold on;
   fill([Z1 Z1 Z2 Z2], [LLR_CutOff, Y1, Y2, LLR_CutOff], clr, 'EdgeColor', clr, 'LineWidth', 1);
end

%% Add vertical lines and titles to each subplot
for nn = 1:nPars
    indPar = Order(nn);
    subplot(nRows, nCols, nn); hold on;
    plot([ParMin(indPar), ParMin(indPar)], [0, LLR_CutOff], 'w:', 'LineWidth', 2);
    set(gca, 'FontName', 'Arial', 'FontSize', 10);

    if ParMin(indPar) < 10
        STR1 = sprintf('%s=%.2g%s', labels{indPar}, ParMin(indPar), Units{indPar});
    else
        STR1 = sprintf('%s=%.0f%s', labels{indPar}, ParMin(indPar), Units{indPar});
    end
    title(STR1, 'FontSize', tfs, 'Color', 'k', 'FontName', 'Arial');
    

    set(gca, 'YDir', 'reverse');
    if nn == indGamma
        set(gca, 'XScale', 'log');
        set(gca, 'XTick', [3 5 10 20 50]);
    end
end
