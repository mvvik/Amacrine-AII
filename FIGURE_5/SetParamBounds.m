function Y = SetParamBounds(X, Rpatch, invFlag)

% "Soft" clipping to upper and lower bounds on parameters:
%
%       X: Parameter array
%  Rpatch: physical domain size (limiting the 1st two parameters)
% invFlag: flag for forward vs. backward "soft-clipping" transform
%
%   Parameters: distances, sensor binding kinetics, and pool fraction:
%      1      2      3        4       5    6        7   
%   Dist-1 Dist-2 SensorKD SensorKp coop gamma pool-fraction

    ParamMin = [  0.002, 0.002,   0.05, 0.001,  1e-8, 0.001, 0.01];
    ParamMax = [ Rpatch, Rpatch,  1000,     1,     2,  1000, 0.99];

    Y = abs( X );

    SoftMin  = @(x, xMin, g) xMin ./ g(xMin ./ x);
    SoftMax  = @(x, xMax, g) xMax .* g(x ./ xMax);

    if invFlag == 0
        Y = SoftMin(Y, ParamMin, @tanh );
        Y = SoftMax(Y, ParamMax, @tanh );
    else
        Y = SoftMax(Y, ParamMax, @atanh );
        Y = SoftMin(Y, ParamMin, @atanh );
    end

    Y = abs(Y);
    Y(isnan(Y)) = ParamMin(isnan(Y));

end

