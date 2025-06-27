% -------------------------------------------------------------------------
%  Likelihood maximization / Parameter profiling for Figures 5B and S4B
% -------------------------------------------------------------------------

nPass   = 4;        	 % Number of profiling passes
nBins   = 36;       	 % Number of slices per parameter (integer times parallel workers)
nTotal  = nPars * nBins; % Number of parameter combinations per single pass

% --- Initialize parallel random streams for each slice
RS = RandStream.create('mrg32k3a', ...
    'NumStreams', nTotal * nPass, ...
    'Seed', 'shuffle', ...
    'CellOutput', true);

% -------------------------------------------------------------------------
%                      Begin profiling passes
% -------------------------------------------------------------------------

for indPass = 1:nPass

    fprintf('\n\n***** Starting pass %d of %d \n\n', indPass, nPass);
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
            Cost0 = 1e6;  cnt = 0;

            while Cost0 > (maxCost + cnt / 2)
                cnt  = cnt + 1;
                X0   = exp(log(Xmin) + (log(Xmax) - log(Xmin)) .* rand(rStream, 1, nPars));
                X0(indPar) = valPar;
                Cost0 = CostFunc(X0);
            end

            fprintf(['Start [Par: %d, Bin: %d, Pass: %d] ERR: %g\n  Pars: ', fStr, '\n\n'], ...
                    indPar, indBin, indPass, Cost0, X0);

            % Profile with fixed parameter:
            X              = removeValue(X0, indPar);                       % Remove current parameter 
            CostFuncSlice  = @(X) CostFunc(insertValue(X, indPar, valPar)); % Cost as function of remaining pars
            [Xopt, Cost]   = fminsearch(CostFuncSlice, X, options);         % Optimize w.r.t. remaining pars
            Y              = insertValue(Xopt, indPar, valPar);             % Insert back current parameter 

            Results(indPar, indBin, :) = [Y, Cost];                         % Store result

            fprintf(['Done [Par: %d, Bin: %d, Pass: %d] %g --> %g\n  Pars: ', fStr, '\n\n'], ...
                    indPar, indBin, indPass, Cost0, Cost, Y);
        end
		fprintf('** Pass %d: Sweep over parameter %d is complete \n\n', indPass, indPar);
    end

    % Reshape and save results; vosualize parameter profile
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


