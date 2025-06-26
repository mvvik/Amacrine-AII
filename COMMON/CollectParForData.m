function ResultsOut = CollectParForData(nPars, DataPrefix)

% Get list of files in the data directory
sls     = ls('../DATA/');
Results = [];

% Loop through all entries and collect matching data
for jj = 1 : size(sls,1)
       if contains(sls(jj,:), DataPrefix ) &&  contains(sls(jj,:), '.mat' ) 
         OutName = ['../DATA/', sls(jj,:)];
         ResultsOut = [];
         load(OutName, 'ResultsOut');
         fprintf('Processing file %s: \n  %d data points\n', sls(jj,:), size(ResultsOut,1));
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
