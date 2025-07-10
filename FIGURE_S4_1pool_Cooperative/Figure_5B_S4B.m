function Figure_5B_S4B()

% -------------------------------------------------------------------------
%   Figure 5B/S4B: Parameter profile plots (max likelihood profile) 
% -------------------------------------------------------------------------

ComputeSetup;
ResultsOut = CollectParForData(nPars, DataPrefix);
if numel(ResultsOut) == 0
    fprintf("No data yet for model %s: run Likelihood_Profile first\n", DataPrefix);
    return;
end
BestError  = 48.3216;  % Smallest error: 2-pool cooperative model with 5 bindings sites

% Parameter labels and units
labels = {'R_X', 'K_D', 'k_{on}', '\beta', '\gamma'};
Units  = {'nm', '\muM', '(\muMms)^{-1}', '', 'ms^{-1}'};

Order    = [1 2 3 4 5];   % Order of parameters for display
nRows    = 2;
nCols    = 3;
nBins    = 30;
clr      = [0.6, 0.0, 0.6];  % Fill color
tfs      = 12;               % Title font size
indGamma = 5;                % Use log scale for parameter #5 = gamma

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocess and filter results

% Rescale first parameter to nm, apply parameter bounds
ResultsOut(:, 1:nPars) = SetParamBounds(ResultsOut(:, 1:nPars), 0);
ResultsOut(:, 1)       = ResultsOut(:, 1) * 1000;
ResultsOut(:, end)     = ResultsOut(:, end) - BestError; 

% Find minimum cost (LLR) and corresponding parameters
[Vmin, Imin] = min(ResultsOut(:, end));
ParMin       = ResultsOut(Imin, 1:nPars);

fprintf(['Minimal error = %g \nRelative to best = %g \n Pars = %g ', fStr, '\n'], Vmin + BestError, Vmin, ParMin);

% Apply cutoff for accepted solutions
MaxError           = Vmin + LLR_CutOff;
validIdx           = ResultsOut(:, end) <= MaxError;
ResultsOut         = ResultsOut(validIdx, :);

fprintf('Filtered set: N = %d\n', size(ResultsOut, 1));

% Reorder parameters for plotting convenience
ResultsOut(:, 1:nPars) = ResultsOut(:, Order);
ParMin                 = ParMin(Order);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Profile likelihood plots for each parameter

figure;

for indPar = 1:nPars
    % Bin parameter values
    Xmin = min(ResultsOut(:, indPar));
    Xmax = max(ResultsOut(:, indPar));
    valPars = zeros(1, nBins);
    valErrs = zeros(1, nBins);

    % Optional log-scaling for parameter 8 (can be skipped here)
    useLog = (indPar == indGamma);
    if useLog
        Xmin = log(Xmin);
        Xmax = log(Xmax);
    end

    for bin = 1:nBins
        % Define bin range
        val0 = Xmin + (bin - 1) / nBins * (Xmax - Xmin);
        val1 = Xmin +  bin      / nBins * (Xmax - Xmin);
        valC = 0.5 * (val0 + val1);

        if useLog
            val0 = exp(val0); val1 = exp(val1); valC = exp(valC);
        end

        valPars(bin) = valC;

        % Extract results within this bin
        inBin = ResultsOut(:, indPar) >= val0 & ResultsOut(:, indPar) < val1;
        if any(inBin)
            valErrs(bin) = min(ResultsOut(inBin, end));
        else
            % Nearest fallback
            [~, closest] = min((ResultsOut(:, indPar) - valC).^2);
            valErrs(bin) = ResultsOut(closest, end);
        end
    end

    % Restore linear range if log-scaling was used
    if useLog
        Xmin = exp(Xmin);
        Xmax = exp(Xmax);
    end

    % Pad ends to close polygon
    valPars = [valPars(1), valPars, valPars(end)];
    valErrs = [Vmin+LLR_CutOff, valErrs, Vmin+LLR_CutOff];

    subplot(nRows, nCols, indPar); hold on;
    fill(valPars, valErrs, clr, 'EdgeColor', clr, 'LineWidth', 1);
    axis([Xmin Xmax Vmin Vmin+LLR_CutOff]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add vertical lines and axis formatting

for nn = 1:nPars
    subplot(nRows, nCols, nn);
    parVal = ParMin(nn);

    % Vertical line at optimal parameter
    plot(parVal * [1 1], [Vmin Vmin+LLR_CutOff], 'w:', 'LineWidth', 2);

    % Annotated title
    if parVal < 10
        strTitle = sprintf('%s = %.2g %s', labels{nn}, parVal, Units{nn});
    else
        strTitle = sprintf('%s = %.0f %s', labels{nn}, parVal, Units{nn});
    end
    title(strTitle, 'FontSize', tfs, 'Color', 'k');

    % Reverse Y-axis to match standard profile-likelihood plot style
    set(gca, 'YDir', 'reverse');

    % Log x-axis scaling for last parameter (assumed gamma)
    if nn == indGamma
        set(gca, 'XScale', 'log');
        set(gca, 'XTick', [3 5 10 20 50]);
    end
    %axis tight;
end
