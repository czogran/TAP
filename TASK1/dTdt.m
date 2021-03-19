function dTdt = dTdt(V,T,delay,Finputs,Tinputs)
    Fh=Finputs(1);
    Fc=Finputs(2);
    Fd=Finputs(3);
    
    Th=Tinputs(1);
    Tc=Tinputs(2);
    Td=Tinputs(3);

    T
    if(delay == 0)
        dTdt=(Fh*Th+Fd*Td-(Fh+Fd)*T)/V;
        return
    else
        dTdt=(Fh*Th+Fc*Tc+Fd*Td-(Fh+Fc+Fd)*T)/V;
    end
end

