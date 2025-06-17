% Save data from "protocol" parfor run, after forced Cntrl-C 

[~, I]   = sort(Results(:,end));
ind0     = find(Results(I,end) > 0, 1);
inds     = I(ind0:end);
indTotal = numel(inds);

ResultsOut = zeros(indTotal, nPars + 1);
ResultsOut(1:indTotal, :) = Results( inds(1:indTotal), : );


save(OutName, 'ResultsOut', '-append');

fprintf('File %s: \n Total data points = %d\n', OutName, indTotal);
