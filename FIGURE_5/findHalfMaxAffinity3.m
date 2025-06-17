function K = findHalfMaxAffinity3( KD, coop )

    % Find half-maximal saturation point for a sensor with 3 binding sites
    K = KD * fsolve( @(x) abs(x).^3 - 3 * coop * abs(x) .* (1 + abs(x)) - 1, 10, optimoptions('fsolve','Display', 'off'));
end

% r^3 / (1 + 3 r beta + 3r^2 beta + r^3) = 1/2
% r^3 - 1 - 3 r beta - 3r^2 beta = 0
% r^3 - 3 r beta (1 + r) = 1