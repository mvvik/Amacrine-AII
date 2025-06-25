function Figure_5B_S4B()

% -------------------------------------------------------------------------
%   Figure 5B/S4B: Parameter profile plots (max likelihood profile) 
% -------------------------------------------------------------------------%

ComputeSetup;
ResultsOut = CollectParForData(nPars);

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

MaxError           = Vmin + LLR_CutOff;
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
for mm = 0 : 6
   ind1 = find(ResultsOut(:,    1) >= 20 + mm*2);
   ind2 = find(ResultsOut(ind1, 1) <  22 + mm*2); 
   inds = ind1(ind2);
   clrG = (0.9 - mm/10);  clrR  = 1 - clrG;  clr = [ clrR, clrG, 0 ];
   sz   = (9 - mm)/8;

    for nn = [1 5] %1 : nPars 
        subplot(nRows, nCols, nn); hold on;
        for ii = 1 : numel(inds)
            XX = ResultsOut(inds(ii),nn); YY = ResultsOut(inds(ii),end);
            plot([XX XX], [YY LLR_CutOff], '-' , 'color', clr, 'linewidth', sz);
        end
    end
end

%% Overlay stratified samples by gamma
for mm = 1 : 6
   ind1 = find(ResultsOut(:,    6) >= (mm + 3)/10);
   ind2 = find(ResultsOut(ind1, 6) <  (mm + 4)/10);
   inds = ind1(ind2);
   clrB = (0.9 - mm/10);  clrR  = 1 - clrB;  clr = [ clrR, 0, clrB ];
   sz   = (9 - mm)/8;

    for nn = [6 8] %1 : nPars 
        subplot(nRows, nCols, nn); hold on;
        for ii = 1 : numel(inds)
            XX = ResultsOut(inds(ii),nn); YY = ResultsOut(inds(ii),end);
            plot([XX XX], [YY LLR_CutOff], '-' , 'color', clr, 'linewidth', sz);
        end
    end
end


%% Add vertical lines and titles to each subplot
for nn = 1:nPars
    indPar = Order(nn);
    subplot(nRows, nCols, nn); hold on;
    plot([ParMin(indPar), ParMin(indPar)], [0, LLR_CutOff], 'w:', 'LineWidth', 2);

    if ParMin(indPar) < 10
        STR1 = sprintf('%s=%.2g%s', labels{indPar}, ParMin(indPar), Units{indPar});
    else
        STR1 = sprintf('%s=%.0f%s', labels{indPar}, ParMin(indPar), Units{indPar});
    end
    title(STR1, 'FontSize', tfs, 'Color', 'k');

    set(gca, 'YDir', 'reverse');
    if nn == indGamma
        set(gca, 'XScale', 'log');
        set(gca, 'XTick', [3 5 10 20 50]);
    end
end
