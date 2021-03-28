function dTdt = dTdtLinearized(V,V0,T,T0,delay,Finputs,Tinputs)
    Fh=Finputs(1);
    Fc=Finputs(2);
    Fd=Finputs(3);
    
    Th=Tinputs(1);
    Tc=Tinputs(2);
    Td=Tinputs(3);

    if(V == 0)
        dTdt = 0;
        return
    end
    if(delay == 0)
        dTdt=(Fh*Th+Fd*Td-(Fh+Fd)*T0)/V0 - (Fh+Fd)*(T-T0)/V0 - ((Fh*Th+Fd*Td-(Fh+Fd)*T0)*(V-V0))/V0^2;
        return
    else
        dTdt=(Fh*Th+Fc*Tc+Fd*Td-(Fh+Fc+Fd)*T0)/V0 - (Fh+Fc+Fd)*(T-T0)/V0 - (Fh*Th+Fc*Tc+Fd*Td-(Fh+Fc+Fd)*T0)/V0^2*(V-V0);
         
    end
end