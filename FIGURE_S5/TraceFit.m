function [Params, X, Y, ERR] = TraceFit(X, Y)
% Fit given Y(X) data to funcion Y(X) = P1 exp(-X/P2) + P3;
% Return the parameters, the Y(X) fit, and the squared log of total error

   ind = find(X > 1, 1);
   X = X(ind:end);
   Y = Y(ind:end);

   Loss = @(P) sum( (log(P(1)*exp(-abs(P(2))*X) ./ X + P(3)) - log(Y)).^2 );

   Params  = fminsearch( Loss, [300 0.02 1], optimset('Display', 'final', 'TolFun', 1e-9));
   Params  = abs( Params );
   
   Y   = Params(1) * exp(-abs(Params(2))*X) ./ X + Params(3);
   ERR = Loss(Params);

end
