function outName = GetNextFileName( prefix )

    sls  = ls('Data/');
    ind  = 1;
    flag = 1;

    while flag
        flag = 0;
        for jj = 1 : size(sls,1)
           %fprintf('%d %d %s %s\n', jj, ind, sls(jj,:), [prefix, num2str(ind), '.mat']) 
           if contains(sls(jj,:), [prefix, num2str(ind), '.mat'] )
             flag = 1;
             ind  = ind + 1;
             break;
           end
        end
    end
    
    %fprintf('%d %d\n', jj, ind) 
    outName = ['Data/', prefix, num2str(ind) ,'.mat'];
    fprintf('Next file: %s\n', outName)
end
