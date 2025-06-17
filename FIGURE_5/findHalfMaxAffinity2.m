function K = findHalfMaxAffinity2( KD, coop )

    % Find half-maximal saturation point for a sensor with 2 binding sites
    K = KD * (coop + sqrt(coop^2 + 1));
end
