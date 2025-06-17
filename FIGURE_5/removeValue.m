function Y = removeValue(X, ind)

    Y      = zeros(1, numel(X) - 1);

    if ind > 1
        Y(1 : (ind-1)) = X(1 : (ind-1));
    end
    if ind < numel(X)
        Y(ind : end) = X((ind + 1) : end);
    end
    
end