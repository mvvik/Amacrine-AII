function EQ = getEquilibriumCoop5( KD, coop )
    
    CaBgr = 0.05; 
    ratio = CaBgr/KD;
    
    S1 = ratio*(5*coop^2);
    S2 = ratio*(2*coop)   * S1;
    S3 = ratio            * S2;
    S4 = ratio/(2*coop)   * S3;
    S5 = ratio/(5*coop^2) * S4;
    EQ = [1, S1, S2, S3, S4, S5];
    EQ = EQ / sum(EQ);
end

