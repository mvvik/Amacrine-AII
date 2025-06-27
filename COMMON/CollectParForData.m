function ResultsOut = CollectParForData(nPars, DataPrefix)
% Collect all Likelihood_Profile output for model defined by nPars and DataPrefix

dataDir = '../DATA/';

% Get list of all .mat files in the directory
files = dir(fullfile(dataDir, '*.mat'));
Results = [];

% Loop through all files
for jj = 1:length(files)
    fileName = files(jj).name;

    % Check for prefix match
    if contains(fileName, DataPrefix)
        OutName = fullfile(dataDir, fileName);
        ResultsOut = [];
        load(OutName, 'ResultsOut');

        fprintf('Processing file %s:\n  %d data points\n', fileName, size(ResultsOut, 1));
        Results = cat(1, Results, ResultsOut);
    end
end
   
% Sort results by the last column and filter out non-positive entries
[~, I]   = sort(Results(:,end));
ind0     = find(Results(I,end) > 0, 1);
inds     = I(ind0:end);
indTotal = numel(inds);

ResultsOut = zeros(indTotal, nPars + 1);

for ii = 1 : indTotal
    ResultsOut(ii, :) = Results( inds(ii), :);
end

% Extract valid sorted results
ResultsOut = abs(ResultsOut); % Ensure values are positive

fprintf('Total data points = %d\n', size(ResultsOut, 1));
