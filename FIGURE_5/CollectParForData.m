% Kappa  = 50;
% BKD    = 1;
% Bkplus = 0.4;
% ICa    = 0.2;               
% Rpatch = 0.250;	

DefineGlobalParameters;

sls     = ls('Data/');
ind     = 1;
flag    = 1;
Results = [];

for jj = 1 : size(sls,1)
       if contains(sls(jj,:), prefixN ) &&  contains(sls(jj,:), '.mat' ) 
         OutName = ['Data/', sls(jj,:)];
         ResultsOut = [];
         load(OutName, 'ResultsOut');
         fprintf('Processing file %s: \n  %d data points\n', sls(jj,:), size(ResultsOut,1));
         %figure; histogram(ResultsOut(:,nPars));
         Results = cat(1, Results, ResultsOut);
       end
end
   
[~, I]   = sort(Results(:,end));
ind0     = find(Results(I,end) > 0, 1);
inds     = I(ind0:end);
indTotal = numel(inds);

ResultsOut = zeros(indTotal, nPars + 1);

for ii = 1 : indTotal
    %ResultsOut(ii, :) = Results( inds(ii), [1 2 3 4 6] );
    ResultsOut(ii, :) = Results( inds(ii), :);
end

ResultsOut = abs(ResultsOut); 

fprintf('Total data points = %d\n', size(ResultsOut, 1));

