function Y = SetParamBounds(X, invFlag)
%
% "Soft" clipping to upper and lower bounds on parameter values:
% ------------------------------------------------------------------
%       X: Parameter array
%  Rpatch: physical domain size (limiting the 1st two parameters)
% invFlag: flag for forward vs. backward "soft-clipping" transform
% ------------------------------------------------------------------
% 1-pool non-cooperative exocytosis sensor parameters:
%     1     2       3     4    
%    R1  SCaKD1  SCakp1 gamma   
    
Rpatch   = 0.25;  % Simulation domain radius
ParamMin = [   0.003,   0.05,  0.001,  0.001];
ParamMax = [  Rpatch,    100,      1,    100];

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
[ii, jj] = find(isnan(Y) | isinf(Y) | abs(imag(Y)) );

for kk = 1 : numel(ii)
    Y(ii(kk),jj(kk)) = 0.5*(ParamMin(jj(kk))+ParamMax(jj(kk)));
end
