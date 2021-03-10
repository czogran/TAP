%% STATIC CHARACTERISTICS %%

clear;
close all;
clc

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
T=33.57;

%VARIABLES
% sample time
Tp=1;
t=0:Tp:2000;

%% static characteristics
%!!!!! IN STATIC CHARACTERISTIC DELAY DOES NOT MATTER
delay=1;

Tvector=[T];

Fc=Fcin;

FhStatic = (Fh-10) : (Fh+10);
FcStatic = (Fc-10) : (Fc+10);

Vstatic = zeros(length(FhStatic), length(FcStatic));
hStatic = zeros(length(FhStatic), length(FcStatic));


for i =  1: length(FhStatic)
    Fh=FhStatic(i);
    for  j=  1: length(FcStatic)
        Fc=FcStatic(j);
        h=81;
        V=[volume(h)];
        Tvector=[T];
        Finputs=[Fh, Fc,Fd];
        for k=2:length(t)
            
            k1= dVdt(h, delay, Finputs);
            k2= dVdt(h+Tp/2*k1,delay,Finputs);
            k3= dVdt(h + Tp/2*k2,delay,Finputs);
            k4= dVdt(h+Tp*k3,delay,Finputs);

            dV=Tp/6*(k1+2*k2+2*k3+k4);
            V(k)=V(k-1)+dV;
    
            h=heightFromVolume(V(k));

%             dVdTdt=Fh*Th+Fc*Tc+Fd*Td-(Fh+Fc+Fd)*Tvector(k-1);
%             Tvector(k)=Tvector(k-1)+dVdTdt/V(k);
        end
        Vstatic(i,j) = V(end);
        hStatic(i,j)=heightFromVolume(V(end));
    end
end

figure
surf(FcStatic, FhStatic,  hStatic)
xlabel('Fh');
ylabel('Fc');
title("Charakterystyka statyczna"+ newline+"wysokości słupa cieczy w zależności od"+newline+"dopływu cieplej i zimnej wody");
hold off