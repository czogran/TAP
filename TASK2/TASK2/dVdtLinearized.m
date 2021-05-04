function dVdt = dVdtLinearized(h,h0, delay, Finputs, vol)
    Fh=Finputs(1);
    Fc=Finputs(2);
    Fd=Finputs(3);

    if(delay == 0)
        if(vol == 1)
            dVdt=Fh+Fd-outputFlowLinearizedVol(volume(h),volume(h0));
        else
            dVdt=Fh+Fd-outputFlowLinearized(h,h0);
        end
    else
        if(vol == 1)
            dVdt=Fh+Fc+Fd-outputFlowLinearizedVol(volume(h),volume(h0));
        else
            dVdt=Fh+Fc+Fd-outputFlowLinearized(h,h0);
        end
    end
end
