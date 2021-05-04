%%--SYMULACJA DZIAŁANIA OBIEKTU--%
%Polecenie: Zasymulować działanie obiektu w Matlabie

clc;
clear;
close all;

%DATA
%F-> dV/dT [cm^3/s]
%T-> temperature

%INPUTS
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

%ADJUSTABLE SIZES
% h;
% Tout;

%CONTROL VALUES
% Fh;
% Fcin;

%VARIABLES
% time
t=0:2000;

%DELAY
delay=120;
delayC=180;

%% filling tank
h=zeros(length(t),1);
V=[0]
Tvector=[0];
Fc=Fcin;
for k=2:length(t)
    if(k<delayC)
        V(k)=V(k-1)+(Fh+Fd-outputFlow(h(k-1)));
        dVdTdt=Fh*Th+Fd*Td-(Fh+Fd)*Tvector(k-1);
    else
        V(k)=V(k-1)+(Fh+Fc+Fd-outputFlow(h(k-1)));
        dVdTdt=Fh*Th+Fc*Tc+Fd*Td-(Fh+Fc+Fd)*Tvector(k-1);
    end
    h(k)=heightFromVolume(V(k));
    
    Tvector(k)=Tvector(k-1)+ dVdTdt/V(k);
end

figure
plot(t,h)
title("Napełnianie zbiornika")
xlabel("t[s]");
ylabel("h[cm]")
hold off

figure
plot(t,Tvector)
title("Temperatura w zbiorniku"+newline+"podczas napełniania zbiornika")
xlabel("t[s]");
ylabel("T[\circC]")
hold off


%% static characteristics
%!!!!! IN STATIC CHARACTERISTIC DELAY DOES NOT MATTER

h=81;
V=[volume(h)]
Tvector=[T];
Fh=Fh;

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
        for k=2:length(t)
            V(k)=V(k-1)+(Fh+Fc+Fd-outputFlow(h));
            h=heightFromVolume(V(k));

            dVdTdt=Fh*Th+Fc*Tc+Fd*Td-(Fh+Fc+Fd)*Tvector(k-1);
            Tvector(k)=Tvector(k-1)+dVdTdt/V(k);
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

function h= heightFromVolume(V)
    %C-> constant
    C=0.7;    
    h=sqrt(V/C);
end

function Tout=outputTemperature(t)
    %delay-> constant;
    delay=120;
    Tout;
end

function Fc=flowC(t)
    %delayC-> constant;
    delayC=180;
    Fc;
end