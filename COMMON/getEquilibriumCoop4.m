function EQ = getEquilibriumCoop4( KD, coop )
    
% Find equilibrium state of a sensor with 4 binding sites

    coop = sqrt( coop );
    CaBgr = 0.05; 
    ratio = CaBgr/KD;
    
    S1 = ratio*(4*coop^3);
    S2 = ratio*(3/2*coop) * S1;
    S3 = ratio/(3/2*coop) * S2;
    S4 = ratio/(4*coop^3) * S3;
    EQ = [1, S1, S2, S3, S4];
    EQ = EQ / sum(EQ);
end

