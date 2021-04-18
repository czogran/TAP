function dVdt = dVdtLinearized(h,h0, delay, Finputs)
    Fh=Finputs(1);
    Fc=Finputs(2);
    Fd=Finputs(3);

    a=7;
    F=a*(sqrt(h0)+(0.5/sqrt(h0))*(h-h0));
    if(delay == 0)
        dVdt=Fh+Fd-F;
    else
        dVdt=Fh+Fc+Fd-F;
    end
end
