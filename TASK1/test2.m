% STEP RESPONSES

%init variables
clc;
clear;
close all;

%DATA
%F-> dV/dT [cm^3/s]
%T-> temperature

%INPUTS
% Fh Fcin Fc Fd
%Hot water
Th=62;
Fh=14;
%Cold water
Tc=23;
Fcin=37;
%Discruption
Td=33;
Fd=12;

%OUTPUT
Tout=0;
F=0;

h=81;
h0=81;
T=33.57;
T0=33.57;

%ADJUSTABLE SIZES
% h;
% Tout;

%CONTROL VALUES
% Fh;
% Fcin-> Fc delayed;

%VARIABLES
% sample time
Tp=1;
t=0:Tp:2000;

%DELAY
delayT=120;
delayC=180;  

% temperature inputs
Tinputs=[Th,Tc,Td];
T0inputs=[Th,Tc,Td];

% inputs in work point
F0inputs=[Fh,Fcin,Fd];

% Vectors of control steps
FcVector =Fcin +[-20,-10,0,10,20]; 
FhVector = Fh +[-14,-10,0,10,20]; 

% auxiliary variables
FcVectorLength = length(FcVector);
FhVectorLength = length(FhVector);

% output vectors of height in tank and linearized height
hVector=ones(length(t),FcVectorLength,FhVectorLength)*h0;
hVectorLinearized=ones(length(t),FcVectorLength,FhVectorLength)*h0;

% volume in work point and output vector of volume in tank
V0=volume(hVector(1));
V=zeros(length(t),FcVectorLength,FhVectorLength);
V(1,:,:)=V0;

% linearized volume in work point and output vector of linearzied volume in tank
VL=zeros(length(t),FcVectorLength,FhVectorLength);
VL(1,:,:)=volume(hVectorLinearized(1));

% temperature in tank and linearized temperature
Tvector=ones(length(t),FcVectorLength,FhVectorLength)*T0;
TvectorL=ones(length(t),FcVectorLength,FhVectorLength)*T0;

% temperature in tank output and linearized temperature
TvectorOutput=ones(length(t),FcVectorLength,FhVectorLength)*T0;
TvectorLOutput=ones(length(t),FcVectorLength,FhVectorLength)*T0;

for j=1:length(FhVector)
    Fh=FhVector(j);
    
    for i=1:length(FcVector)    
    
    % rungy-kutta
    for k=2:length(t)
        delay=1;

        if(k<delayC)
            Fc=Fcin;
        else
           Fc=FcVector(i);
        end 
        Finputs=[Fh,Fc,Fd];
         
        h=hVector(k-1,i,j);
        hL=hVectorLinearized(k-1,i,j);
        
        T=Tvector(k-1,i,j);
        TL=TvectorL(k-1,i,j);
        
        kV1= dVdt(h, delay, Finputs);
        kV2= dVdt(heightFromVolume(V(k-1,i,j) + Tp/2*kV1),delay,Finputs);
        kV3= dVdt(heightFromVolume(V(k-1,i,j) + Tp/2*kV2),delay,Finputs);
        kV4= dVdt(heightFromVolume(V(k-1,i,j) + Tp*kV3),delay,Finputs);

        dV=Tp/6*(kV1+2*kV2+2*kV3+kV4);
        V(k,i,j)=V(k-1,i,j)+dV;
        hVector(k,i,j)=heightFromVolume(V(k,i,j));

        kV1L= dVdtLinearized(hL,h0, delay, Finputs);
        kV2L= dVdtLinearized(heightFromVolume(VL(k-1,i,j) + Tp/2*kV1L),h0,delay,Finputs);
        kV3L= dVdtLinearized(heightFromVolume(VL(k-1,i,j) + Tp/2*kV2L),h0,delay,Finputs);
        kV4L= dVdtLinearized(heightFromVolume(VL(k-1,i,j) + Tp*kV3L),h0,delay,Finputs);

        dVL=Tp/6*(kV1L+2*kV2L+2*kV3L+kV4L);
        VL(k,i,j)=VL(k-1,i,j)+dVL;
        hVectorLinearized(k,i,j)=heightFromVolume(VL(k,i,j));
        
                
        kT1= dTdt(V(k-1,i,j),T,delay,Finputs,Tinputs);
        kT2= dTdt(V(k-1,i,j) + Tp*V(k-1,i,j)/2,T + Tp/2*kT1,delay,Finputs,Tinputs);
        kT3= dTdt(V(k-1,i,j) + Tp*V(k-1,i,j)/2,T + Tp/2*kT2,delay,Finputs,Tinputs);
        kT4= dTdt(V(k-1,i,j) + Tp*V(k-1,i,j),T + Tp*kT3,delay,Finputs,Tinputs);       

        dT=Tp/6*(kT1+2*kT2+2*kT3+kT4);
        Tvector(k,i,j)=Tvector(k-1,i,j)+ dT;
        

        kT1L= dTdtLinearized(VL(k-1,i,j),V0,TL,T0,delay,Finputs,F0inputs,Tinputs,T0inputs);
        kT2L= dTdtLinearized(VL(k-1,i,j) + Tp/2*kV1L,V0,TL + Tp/2*kT1L,T0,delay,Finputs,F0inputs,Tinputs,T0inputs);
        kT3L= dTdtLinearized(VL(k-1,i,j) + Tp/2*kV2L,V0,TL + Tp/2*kT2L,T0,delay,Finputs,F0inputs,Tinputs,T0inputs);
        kT4L= dTdtLinearized(VL(k-1,i,j) + Tp*kV3L,V0,TL + Tp*kT3L,T0,delay,Finputs,F0inputs,Tinputs,T0inputs);
        
        dTLinearized=Tp/6*(kT1L+2*kT2L+2*kT3L+kT4L);
        TvectorL(k,i,j)=TvectorL(k-1,i,j)+dTLinearized;
        
        if(k<=delayT)
            TvectorOutput(k,i,j)=T0;
            TvectorLOutput(k,i,j)=T0;
        else
            TvectorOutput(k,i,j)=Tvector(k-delayT,i,j);
            TvectorLOutput(k,i,j)=TvectorL(k-delayT,i,j);
        end
    end
end
end

