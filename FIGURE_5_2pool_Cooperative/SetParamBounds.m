function Y = SetParamBounds(X, invFlag)
%
% "Soft" clipping to upper and lower bounds on parameter values:
% ----------------------------------------------------------------
%       X: Parameter array
% invFlag: flag for forward vs. backward "soft-clipping" transform
% ----------------------------------------------------------------
% 2-pool cooperative exocytosis sensor parameters:
%      1      2      3        4       5    6        7   
%   Dist-1 Dist-2 SensorKD SensorKp coop gamma pool-fraction

Rpatch = 0.25;      % Simulation domain radius
Y      = abs( X );  % Discard negative signs on parameters

ParamMin = [  0.003, 0.003,   0.05,  0.001,  1e-5, 0.001, 0.01];
ParamMax = [ Rpatch, Rpatch,   100,      1,     2,   100, 0.99];

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
[ii, jj] = find(isnan(Y) | isinf(Y) | abs(imag(Y)) );

for kk = 1 : numel(ii)
    Y(ii(kk),jj(kk)) = 0.5*(ParamMin(jj(kk))+ParamMax(jj(kk)));
end


