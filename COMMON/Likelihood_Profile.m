% -------------------------------------------------------------------------
%  Likelihood maximization / Parameter profiling for Figures 5B and S4B
% -------------------------------------------------------------------------
%  Note: algorithm could be optimized by sharing information between passes
% -------------------------------------------------------------------------

if contains(pwd,'COMMON')
    fprintf('Do not execute scripts in this folder (switch to a Figure folder)\n');
    return;
end

% -------------------------------------------------------------------------

nPass   = 5;        	 % Number of profiling passes
nBins   = 36;       	 % Number of slices per parameter (integer times parallel workers)
nTotal  = nPars * nBins; % Number of parameter combinations per single pass

fprintf('Likelihood profile computation:\n');
fprintf('Model %s with %d calcium binding sites: %d parameters\n', DataPrefix, nCaSites, nPars);

% --- Initialize parallel random streams for each slice
RS = RandStream.create('mrg32k3a', ...
    'NumStreams', nTotal * nPass, ...
    'Seed', 'shuffle', ...
    'CellOutput', true);

% -------------------------------------------------------------------------
%                      Begin profiling passes
% -------------------------------------------------------------------------

for indPass = 1:nPass

    fprintf('\n***** Starting pass %d of %d \n\n', indPass, nPass);
    Results = zeros(nPars, nBins, nPars + 1);   % Allocate results tensor

    % Optimization settings
    options = optimset('Display', 'off', 'TolFun', 3e-3, 'TolX', 1e-4, 'MaxFunEvals', 500);

    % Create output file for this pass
    OutName = ['../DATA/', GetNextFileName(DataPrefix)];
    save(OutName, 'nPars');

    % Loop over parameters to fix each in turn
    for indPar = 1:nPars

        % Loop over bins for the fixed parameter
        parfor indBin = 1:nBins

            % Compute the fixed parameter value for this bin
            valPar = Xmin(indPar) + (indBin - 1) / (nBins - 1) * (Xmax(indPar) - Xmin(indPar));

            fprintf("Par: %d, Bin: %d, Pass: %d, val: %g \n", indPar, indBin, indPass, valPar);

            % Use a unique random stream
            rStream = RS{indBin * indPar * indPass};

            % Randomly sample an initial guess until it's good enough
            cnt = 0;  MinCost = 1e80; MinPar = zeros(nPars, 1);

            while MinCost > (maxCost + cnt/2)  % Becomes less greedy with each iteration
                cnt  = cnt + 1;
                X0   = exp(log(Xmin) + (log(Xmax) - log(Xmin)) .* rand(rStream, 1, nPars));
                X0(indPar) = valPar;
                Cost0 = CostFunc(X0);
                if Cost0 < MinCost
                    MinCost = Cost0;
                    MinPar  = X0;
                end
            end

            fprintf(['Start [Par: %d, Bin: %d, Pass: %d] ERR: %g\n  Pars: ', fStr, '\n\n'], ...
                    indPar, indBin, indPass, MinCost, MinPar);

            % Profile with fixed parameter:
            CostFuncSlice  = @(X) CostFunc(insertValue(X, indPar, valPar)); % Cost as function of remaining pars
            X0             = removeValue(MinPar, indPar);                   % Remove current parameter      
            [Xopt, Cost]   = fminsearch(CostFuncSlice, X0, options);        % Optimize w.r.t. remaining pars
            Y              = insertValue(Xopt, indPar, valPar);             % Insert back current parameter 

            Results(indPar, indBin, :) = [Y, Cost];                         % Store result

            fprintf(['Done [Par: %d, Bin: %d, Pass: %d] %g --> %g\n  Pars: ', fStr, '\n\n'], ...
                    indPar, indBin, indPass, MinCost, Cost, Y);
        end
		fprintf('** Pass %d: Sweep over parameter %d done (min = %g)\n\n', ...
		         indPass, indPar, min(Results(indPar,:,end)));
    end

    % Reshape and save results; visualize parameter profile
    ResultsOut = reshape(Results, nTotal, nPars + 1);
    SaveParForData(ResultsOut, nPars, OutName);
    Figure_5B_S4B;
    drawnow;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Y = insertValue(X, ind, val)
% Insert a value into vector X at position ind, shifting the rest
% Example: insertValue([1 2 3], 2, 99) -> [1 99 2 3]
    Y = [X(1:ind-1), val, X(ind:end)];
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Y = removeValue(X, ind)
%REMOVEVALUE Remove value at given index from vector X.
    Y = [X(1:ind-1), X(ind+1:end)];
end


