function [rArray, CaArray, tArray] = ReadDataVsTime(fName)

    %fprintf('Attempting to read %s \n', fName);
    f      = fopen(fName, 'rb');      % Read the dump file
    G      = fread(f, 1,  'int');     % Read geometry (lowest 2 bits = # of dimensions)
    N      = fread(f, 1,  'int');     % Read number of r-nodes in binary plot
    rArray = fread(f, N,  'double');  % r-axis coordinates
    
    state   = 0;
    CaArray = [];
    tArray  = [];
    tOld    = -1;
    
    while ~state
        tNew = fread(f, 1, 'double');  % Read Time
        Ca   = fread(f, N, 'double');  % Read the data
        
        if feof(f) 
            fclose(f);
            return;  
        end
    
        if tNew > tOld
            tArray  = [ tArray, tNew ];
            CaArray = [ CaArray, Ca  ];
            tOld = tNew;
        end
        
        state = feof(f);
    end
    
    fclose(f);
end