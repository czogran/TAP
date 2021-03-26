%%--SYMULACJA DZIAŁANIA OBIEKTU--%
%Polecenie: Zasymulować działanie obiektu w Matlabie

%init data and inports
init

%% tank filling
hVector=zeros(length(t),1);

V=[];
V = zeros(1, 2001);
V(1) = 0;
TVector=zeros(1, 2001);
TVector(1) = 0;
ToutputVector=zeros(1, 2001);

Fc=Fcin;
Finputs=[Fh,Fc,Fd];
Tinputs=[Th,Tc,Td];

% Euler "normal"
for k=2:length(t)
    h=hVector(k-1);
    T=TVector(k-1);
    
    if(k<delayC)
      delayFc=0;  
    else
       delayFc=1;
    end
    
   
    dV=Tp*dVdt(h,delayFc, Finputs);
    V(k)=V(k-1)+dV;
    hVector(k)=heightFromVolume(V(k));
    
%     TODO-czy dobre indeksy??????
    dT=dTdt(V(k-1),T,delayFc,Finputs,Tinputs);
    TVector(k)=TVector(k-1)+dT;
    
    
    if(k<=delayT)
        ToutputVector(k)=0;
    else
        ToutputVector(k)=TVector(k-delayT);
    end
end

eNormalT = TVector;
eNormalTout = ToutputVector;

figure
plot(t,hVector)
title("Napełnianie zbiornika" + newline + "symulacja metodą Eulera")
xlabel("t[s]");
ylabel("h[cm]")
hold off

figure
plot(t,TVector)
hold on
plot(t,ToutputVector)
title("Napełnianie zbiornika" + newline+"tempteratura"+newline + "symulacja metodą Eulera")
xlabel("t[s]");
ylabel("T[\circC]")
legend("temperatura w zbiorniku", "temperatura na wyjsciu", 'Location','best')
hold off

eNormal=hVector;

% Euler modified
for k=2:length(t)
    h=hVector(k-1);
    T=TVector(k-1);
    
    if(k<delayC)
      delay=0;  
    else
       delay=1;
    end
    
   pomV=V(k-1)+0.5*Tp*dVdt(h,delay, Finputs);
   pomh=heightFromVolume(pomV);
   
   dV=Tp*dVdt(pomh,delay, Finputs);
   V(k)=V(k-1)+dV;
   
   pomT=T + 0.5*Tp*dTdt(V(k-1),T,delayFc,Finputs,Tinputs);
   dT=Tp*dTdt(pomV,pomT,delayFc,Finputs,Tinputs);
   TVector(k)=TVector(k-1)+dT;
    
    
   if(k<=delayT)
       ToutputVector(k)=0;
   else
       ToutputVector(k)=TVector(k-delayT);
   end
   
   
   hVector(k)=heightFromVolume(V(k));
   
   
   
%     Tvector(k)=Tvector(k-1)+ dVdTdt/V(k);
end

eModifiedT = TVector;
eModifiedTout = ToutputVector;

figure
plot(t,hVector)
title("Napełnianie zbiornika" + newline + "symulacja zmodyfikowaną metodą Eulera")
xlabel("t[s]");
ylabel("h[cm]")
hold off

eModified = hVector;

% figure
% plot(t,Tvector)
% title("Temperatura w zbiorniku"+newline+"podczas napełniania zbiornika")
% xlabel("t[s]");
% ylabel("T[\circC]")
% hold off
h0 = 80;
V0 = volume(h0);
T0 = 30;
hVector=ones(1, length(t))*0;
hVectorLinearized=ones(1, length(t))*0;

V=volume(hVector);
VLinearized=volume(hVectorLinearized);

TVector=ones(1, length(t))*0;
TVectorLinearized=ones(1, length(t))*0;

doVol = 0;

% rungy-kutta
for k=2:length(t)
    h=hVector(k-1);
    T=TVector(k-1);
    TL=TVectorLinearized(k-1);
    hL=hVectorLinearized(k-1);
    
    if(k<delayC)
      delay=0;  
    else
       delay=1;
    end
    kV1= dVdt(h, delay, Finputs);
    kV2= dVdt(heightFromVolume(V(k-1) + Tp/2*kV1),delay,Finputs);
    kV3= dVdt(heightFromVolume(V(k-1) + Tp/2*kV2),delay,Finputs);
    kV4= dVdt(heightFromVolume(V(k-1) + Tp*kV3),delay,Finputs);
    
    kV1L= dVdtLinearized(hL,h0, delay, Finputs, doVol);
    kV2L= dVdtLinearized(heightFromVolume(VLinearized(k-1)+Tp/2*kV1L),h0,delay,Finputs, doVol);
    kV3L= dVdtLinearized(heightFromVolume(VLinearized(k-1)+Tp/2*kV2L),h0,delay,Finputs, doVol);
    kV4L= dVdtLinearized(heightFromVolume(VLinearized(k-1)+Tp*kV3L),h0,delay,Finputs, doVol);
    
    kT1= dTdt(V(k-1),T,delay,Finputs,Tinputs);
    kT2= dTdt(V(k-1) + Tp/2*kV1,T + Tp/2*kT1,delay,Finputs,Tinputs);
    kT3= dTdt(V(k-1) + Tp/2*kV2,T + Tp/2*kT2,delay,Finputs,Tinputs);
    kT4= dTdt(V(k-1) + Tp*kV3,T + Tp*kT3,delay,Finputs,Tinputs);
    
    kT1L= dTdtLinearized(VLinearized(k-1),V0,TL,T0,delay,Finputs,Tinputs);
    kT2L= dTdtLinearized(VLinearized(k-1) + Tp/2*kV1L,V0,TL + Tp/2*kT1L,T0,delay,Finputs,Tinputs);
    kT3L= dTdtLinearized(VLinearized(k-1) + Tp/2*kV2L,V0,TL + Tp/2*kT2L,T0,delay,Finputs,Tinputs);
    kT4L= dTdtLinearized(VLinearized(k-1) + Tp*kV3L,V0,TL + Tp*kT3L,T0,delay,Finputs,Tinputs);
    
    dV=Tp/6*(kV1+2*kV2+2*kV3+kV4);
    V(k)=V(k-1)+dV;
    
    dVLinearized=Tp/6*(kV1L+2*kV2L+2*kV3L+kV4L);
    VLinearized(k)=VLinearized(k-1)+dVLinearized;
    
    dT=Tp/6*(kT1+2*kT2+2*kT3+kT4);
    dTLinearized=Tp/6*(kT1L+2*kT2L+2*kT3L+kT4L);
    
    hVector(k)=heightFromVolume(V(k));
    hVectorLinearized(k)=heightFromVolume(VLinearized(k));
    
    %dT=dTdt(V(k),T,delayFc,Finputs,Tinputs);
    TVector(k)=TVector(k-1)+dT;
    TVectorLinearized(k)=TVectorLinearized(k-1)+dTLinearized;
    
    
    if(k<=delayT)
        ToutputVector(k)=0;
    else
        ToutputVector(k)=TVector(k-delayT);
    end
%     Tvector(k)=Tvector(k-1)+ dVdTdt/V(k);



end

rkT = TVector;
rkTout = ToutputVector;


rkTLinearized = TVectorLinearized;

figure
plot(t,hVector)
title("Napełnianie zbiornika" + newline + "symulacja metodą Rungego Kutty")
hold on
plot(t, hVectorLinearized)
xlabel("t[s]");
ylabel("h[cm]")
hold off

rk = hVector;


figure
plot(t,eNormal, 'r')
hold on
plot(t,eModified,'g')
hold on
plot(t, rk, 'b')
title("Napełnianie zbiornika" + newline + "symulacja metodą trzech różnych metod")
legend("zwykła metoda Eulera", "zmodyfikowana metoda Eulera", "metoda Rungego Kutty", 'Location', 'best');
xlabel("t[s]");
ylabel("h[cm]")
hold off

figure
plot(t,rkT)
hold on
plot(t,rkTLinearized)
title("Napełnianie zbiornika" + newline+"tempteratura"+newline + "symulacja metodą RK4")
xlabel("t[s]");
ylabel("T[\circC]")
legend("temperatura w zbiorniku", "temperatura zlinearyzowana", 'Location','best')
hold off


dlmwrite('DefModel_NorE.txt', [t' eNormal], 'delimiter', ' ')
dlmwrite('DefModel_ModE.txt', [t' eModified], 'delimiter', ' ')
dlmwrite('DefModel_RK4.txt', [t rk], 'delimiter', ' ')

dlmwrite('DefTModel_NorE.txt', [t' eNormalT'], 'delimiter', ' ')
dlmwrite('DefTModel_ModE.txt', [t' eModifiedT'], 'delimiter', ' ')
dlmwrite('DefTModel_RK4.txt', [t' rkT'], 'delimiter', ' ')

dlmwrite('DefToutModel_NorE.txt', [t' eNormalTout'], 'delimiter', ' ')
dlmwrite('DefToutModel_ModE.txt', [t' eModifiedTout'], 'delimiter', ' ')
dlmwrite('DefToutModel_RK4.txt', [t' rkTout'], 'delimiter', ' ')

name="method-comparison-flow";
fileName="object-simulation";
overLeafFilePath="img/object-simulation/";
caption="Symulacja napełniania zbironika przy użyciu trzech różnych metod obliczeń numerycznych";
label="fig:method-comparison-flow";

saveFigure(gcf, 'img/object-simulation/',name,fileName,overLeafFilePath,caption,label);