function Y = insertValue(X, ind, val)

    Y      = zeros(1, numel(X) + 1);
    Y(ind) = val;

    if ind > 1
        Y(1 : (ind-1)) = X(1 : (ind-1));
    end
    if ind <= numel(X) 
        Y((ind + 1) : end) = X(ind : end);
    end
    
end