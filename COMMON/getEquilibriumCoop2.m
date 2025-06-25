function EQ = getEquilibriumCoop2( KD, coop )
    
% Find equilibrium state of a sensor with 2 binding sites
    
    coop = sqrt( coop );
    CaBgr = 0.05; 
    ratio = CaBgr/KD;
    
    S1 = 2*coop*ratio;    
    S2 = ratio^2;
    
    EQ = [1, S1, S2] ./ (1 + S1 + S2);
end

