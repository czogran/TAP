function dVdt = dVdt(h, delay, Finputs)
    Fh=Finputs(1);
    Fc=Finputs(2);
    Fd=Finputs(3);

    if(delay == 0)
        dVdt=Fh+Fd-outputFlow(h);
    else
        dVdt=Fh+Fc+Fd-outputFlow(h);
    end
end
