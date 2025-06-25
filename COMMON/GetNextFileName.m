function outName = GetNextFileName( prefix )
% Create another file name with given prefix, by adding an integer index

    sls  = ls('../DATA/');
    ind  = 1;
    flag = 1;

    while flag
        flag = 0;
        for jj = 1 : size(sls,1)
           if contains(sls(jj,:), [prefix, num2str(ind), '.mat'] )
             flag = 1;
             ind  = ind + 1;
             break;
           end
        end
    end
    
    outName = ['../DATA/', prefix, num2str(ind) ,'.mat'];
    fprintf('Next file: %s\n', outName)
end
