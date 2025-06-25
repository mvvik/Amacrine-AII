function EQ = Get_Equilibrium( KD1, KD2 )

    CaBgr = 0.05;
    
    ratio1 = CaBgr / KD1;
    ratio2 = CaBgr / KD2 * ratio1;
    
    S0 = 1 / (1 + ratio1 + ratio2);
    S1 = ratio1 * S0;
    S2 = ratio2 * S0;

    EQ = [S0, S1, S2];
end

