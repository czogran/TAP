function dVdt = dVdtLinearized(h,h0, delay, Finputs)
    Fh=Finputs(1);
    Fc=Finputs(2);
    Fd=Finputs(3);

    if(delay == 0)
        dVdt=Fh+Fd-outputFlowLinearized(h,h0);
    else
        dVdt=Fh+Fc+Fd-outputFlowLinearized(h,h0);
    end
end
