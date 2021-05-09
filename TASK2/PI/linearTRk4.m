kT1= dTdtLinearized(V,V0,T,T0,delay,Finputs,F0inputs,Tinputs,T0inputs);
kT2= dTdtLinearized(V + Ts/2*kV1,V0,TL + Tp/2*kT1,T0,delay,Finputs,F0inputs,Tinputs,T0inputs);
kT3= dTdtLinearized(V + Ts/2*kV2,V0,TL + Tp/2*kT2,T0,delay,Finputs,F0inputs,Tinputs,T0inputs);
kT4= dTdtLinearized(V + Ts*kV3,V0,TL + Tp*kT3L,T0,delay,Finputs,F0inputs,Tinputs,T0inputs);
