%%--SYMULACJA DZIAŁANIA OBIEKTU--%
%Polecenie: Zasymulować działanie obiektu w Matlabie

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
T=33.57;

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
delay=120;
delayC=180;

%% tank filling
hVector=zeros(length(t),1);
V=[0];
Tvector=[0];
Fc=Fcin;
Finputs=[Fh,Fc,Fd];

% Euler zwykły
for k=2:length(t)
    h=hVector(k-1);

    if(k<delayC)
      delay=0;  
    else
       delay=1;
    end
    
   
    dV=Tp*dVdt(h,delay, Finputs);
    V(k)=V(k-1)+dV;
    hVector(k)=heightFromVolume(V(k));
    
%     Tvector(k)=Tvector(k-1)+ dVdTdt/V(k);
end

figure
plot(t,hVector)
title("Napełnianie zbiornika" + newline + "symulacja metodą Eulera")
xlabel("t[s]");
ylabel("h[cm]")
hold off

e=hVector;

% Euler zmodyfikowany
for k=2:length(t)
    h=hVector(k-1);

    if(k<delayC)
      delay=0;  
    else
       delay=1;
    end
    
   pomV=V(k-1)+0.5*Tp*dVdt(h,delay, Finputs);
   pomh=heightFromVolume(pomV);
   
   dV=Tp*dVdt(pomh,delay, Finputs);
   V(k)=V(k-1)+dV;
     hVector(k)=heightFromVolume(V(k));
    
%     Tvector(k)=Tvector(k-1)+ dVdTdt/V(k);
end

figure
plot(t,hVector)
title("Napełnianie zbiornika" + newline + "symulacja zmodyfikowaną metodą Eulera")
xlabel("t[s]");
ylabel("h[cm]")
hold off

% figure
% plot(t,Tvector)
% title("Temperatura w zbiorniku"+newline+"podczas napełniania zbiornika")
% xlabel("t[s]");
% ylabel("T[\circC]")
% hold off

% rungy-kutta
for k=2:length(t)
    h=hVector(k-1);

    if(k<delayC)
      delay=0;  
    else
       delay=1;
    end
    k1= dVdt(h, delay, Finputs);
    k2= dVdt(h+Tp/2*k1,delay,Finputs);
    k3= dVdt(h + Tp/2*k2,delay,Finputs);
    k4= dVdt(h+Tp*k3,delay,Finputs);
    
    dV=Tp/6*(k1+2*k2+2*k3+k4);
    V(k)=V(k-1)+dV;
    
    hVector(k)=heightFromVolume(V(k));
    
%     Tvector(k)=Tvector(k-1)+ dVdTdt/V(k);
end

figure
plot(t,hVector)
title("Napełnianie zbiornika" + newline + "symulacja metodą Rungego Kutty")
xlabel("t[s]");
ylabel("h[cm]")
hold off


% %% static characteristics
% %!!!!! IN STATIC CHARACTERISTIC DELAY DOES NOT MATTER
% 
% h=81;
% V=[volume(h)]
% Tvector=[T];
% Fh=Fh;
% 
% FhStatic = (Fh-10) : (Fh+10);
% FcStatic = (Fc-10) : (Fc+10);
% 
% Vstatic = zeros(length(FhStatic), length(FcStatic));
% hStatic = zeros(length(FhStatic), length(FcStatic));
% 
% for i =  1: length(FhStatic)
%     Fh=FhStatic(i);
%     for  j=  1: length(FcStatic)
%         Fc=FcStatic(j);
%         h=81;
%         V=[volume(h)];
%         Tvector=[T];
%         for k=2:length(t)
%             V(k)=V(k-1)+(Fh+Fc+Fd-outputFlow(h));
%             h=heightFromVolume(V(k));
% 
%             dVdTdt=Fh*Th+Fc*Tc+Fd*Td-(Fh+Fc+Fd)*Tvector(k-1);
%             Tvector(k)=Tvector(k-1)+dVdTdt/V(k);
%         end
%         Vstatic(i,j) = V(end);
%         hStatic(i,j)=heightFromVolume(V(end));
%     end
% end
% 
% figure
% surf(FcStatic, FhStatic,  hStatic)
% xlabel('Fh');
% ylabel('Fc');
% title("Charakterystyka statyczna"+ newline+"wysokości słupa cieczy w zależności od"+newline+"dopływu cieplej i zimnej wody");
% hold off



