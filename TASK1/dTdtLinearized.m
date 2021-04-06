function dTdt = dTdtLinearized(V,V0,T,T0,delay,Finputs,F0inputs,Tinputs, T0inputs)
    Fh0=F0inputs(1);
    Fc0=F0inputs(2);
    Fd0=F0inputs(3);

    Fh=Finputs(1);
    Fc=Finputs(2);
    Fd=Finputs(3);
    
    Th=Tinputs(1);
    Tc=Tinputs(2);
    Td=Tinputs(3);
    
    Th0=T0inputs(1);
    Tc0=T0inputs(2);
    Td0=T0inputs(3);

    if(V == 0)
        dTdt = 0;
        return
    end
    
    

    
    if(delay == 0)
        %dTdt=(Fh0*Th+Fd0*Td-(Fh0+Fd0)*T0)/V0-...
        %(Fh0+Fd0)/V0*(T-T0)+...
        %((Fh0+Fd0)*T0-(Fh0*Th+Fd0*Td))/V0^2*(V-V0)+...
        %(Th-T0)/V0*(Fh-Fh0)+...
        %(Td-T0)/V0*(Fd-Fd0);
%         dTdt=(Fh*Th+Fd*Td-(Fh+Fd)*T0)/V0 - (Fh+Fd)*(T-T0)/V0 - ((Fh*Th+Fd*Td-(Fh+Fd)*T0)*(V-V0))/V0^2;
        dTdt = (Th0 - T0)/V0*Fh + (Td0 - T0)/V0*Fd + Fh0/V0*Th + Fd0/V0*Td + (Fh0*(T0 - Th0) + Fd0*(T0 - Td0))/(V0^2)*V - (Fh0 + Fd0)/V0*T;
        return
    else
        %dTdt=(Fh0*Th+Fd0*Td+Fc0*Tc-(Fh0+Fd0+Fc0)*T0)/V0-...
        %(Fh0+Fc0+Fd0)/V0*(T-T0)+...
        %((Fh0+Fc0+Fd0)*T0-(Fh0*Th+Fc0*Tc+Fd0*Td))/V0^2*(V-V0)+...
        %(Th-T0)/V0*(Fh-Fh0)+...
        %(Tc-T0)/V0*(Fc-Fc0)+...
        %(Td-T0)/V0*(Fd-Fd0);
        dTdt = (Th0 - T0)/V0*Fh + (Tc0 - T0)/V0*Fc + (Td0 - T0)/V0*Fd + Fh0/V0*Th + Fc0/V0*Tc + Fd0/V0*Td + (Fh0*(T0 - Th0) + Fc0*(T0 - Tc0) + Fd0*(T0 - Td0))/(V0^2)*V - (Fh0 + Fc0 + Fd0)/V0*T;
        
%         dTdt=(Fh*Th+Fc*Tc+Fd*Td-(Fh+Fc+Fd)*T0)/V0 - (Fh+Fc+Fd)*(T-T0)/V0 - (Fh*Th+Fc*Tc+Fd*Td-(Fh+Fc+Fd)*T0)/V0^2*(V-V0);
         
    end
end