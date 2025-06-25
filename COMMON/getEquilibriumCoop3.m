function EQ = getEquilibriumCoop3( KD, coop )
    
% Find equilibrium state of a sensor with 3 binding sites

    CaBgr = 0.05; 
    ratio = CaBgr/KD;
    
    S1 = ratio*3*coop;
    S2 = ratio          * S1;
    S3 = ratio/(3*coop) * S2;
    EQ = [1, S1, S2, S3];
    EQ = EQ / sum(EQ);
end

