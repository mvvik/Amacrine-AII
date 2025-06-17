function K = findHalfMaxAffinity4( KD, coop )

    % Find half-maximal saturation point for a sensor with 4 binding sites
    K = KD * fsolve( @(x) abs(x).^4 - 2 * coop^3 * abs(x) .* (2 + 3*coop*abs(x) + 2*abs(x).^2) - 1, 10, optimoptions('fsolve','Display', 'off'));
end
