function dVdt = dVdt(h, delay, Finputs)
    Fh=Finputs(1);
    Fc=Finputs(2);
    Fd=Finputs(3);

    a=7;
    F=a*sqrt(h);
    if(delay == 0)
        dVdt=Fh+Fd-F;
    else
        dVdt=Fh+Fc+Fd-F;
    end
end
