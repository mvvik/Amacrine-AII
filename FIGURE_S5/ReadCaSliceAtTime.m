function [Ca, rArray]  = ReadCaSliceAtTime( fName, newTime )

    [rArray, CaGrid, tArray] = ReadDataVsTime(fName);
    [TT, RR] = meshgrid(tArray,  rArray );
    [tt, rr] = meshgrid(newTime, rArray );
    Ca = interp2(TT, RR, CaGrid, tt, rr, 'spline');
    
end