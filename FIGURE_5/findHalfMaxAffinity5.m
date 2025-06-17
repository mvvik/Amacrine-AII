function K = findHalfMaxAffinity5( KD, coop )

    % Find half-maximal saturation point for a sensor with 5 binding sites
     K = KD * fsolve( @(x) abs(x).^5 - 5 * coop^2 * abs(x) .* (1 + 2*coop*abs(x).*(1 + abs(x)) + abs(x).^3) - 1, 10, optimoptions('fsolve','Display', 'off'));
end
