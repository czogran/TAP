kT1= dTdt(V,T,delay,Finputs,Tinputs);
kT2= dTdt(V + Ts*V/2,T + Ts/2*kT1,delay,Finputs,Tinputs);
kT3= dTdt(V + Ts*V/2,T + Ts/2*kT2,delay,Finputs,Tinputs);
kT4= dTdt(V + Ts*V,T + Ts*kT3,delay,Finputs,Tinputs);       
