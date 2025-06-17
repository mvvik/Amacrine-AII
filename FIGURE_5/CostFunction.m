function [cost, Bnd] = CostFunction(GetModel, GetEquil, X, timeArray1, timeArray2, timeArray3, CaGrid1, CaGrid2, CaGrid3, DataIn, ErrorIn, rArray, DT)

    DefineLocalParameters;

    X = SetParamBounds( X, Rpatch, 0 );
   
    EQ  = GetEquil( X(3), X(5) );
    Bnd = EQ(end);

    odeSolve = @ode15s;
	ACC      = 1e-5;

    Out1 = GetModel( EQ, CaGrid1, timeArray1, rArray, X, DT, odeSolve, ACC );
    Out2 = GetModel( EQ, CaGrid2, timeArray2, rArray, X, DT, odeSolve, ACC );
    Out3 = GetModel( EQ, CaGrid3, timeArray3, rArray, X, DT, odeSolve, ACC );
    
    DataOut = [Out1, Out2, Out3];
    scale = sum(DataIn .* DataOut .* ErrorIn) / sum(DataOut.* DataOut .* ErrorIn);
    cost  = sum( (DataIn - scale * DataOut).^2 .* ErrorIn) * (1 + 20 * Bnd^2);
end

