%%--SYMULACJA DZIAŁANIA OBIEKTU--%
%Polecenie: Zasymulować działanie obiektu w Matlabie

%init data and inports
init

%% tank filling
hVector=zeros(length(t),1);

V=[0];
Tvector=[0];
ToutputVector=[0];

Fc=Fcin;
Finputs=[Fh,Fc,Fd];
Tinputs=[Th,Tc,Td];

% Euler "normal"
for k=2:length(t)
    h=hVector(k-1);
    T=Tvector(k-1);
    
    if(k<delayC)
      delayFc=0;  
    else
       delayFc=1;
    end
    
    dV=Tp*dVdt(h,delayFc, Finputs);
    V(k)=V(k-1)+dV;
    hVector(k)=heightFromVolume(V(k));
    
%     TODO-czy dobre indeksy??????
    dT=dTdt(V(k),T,delayFc,Finputs,Tinputs);
    Tvector(k)=Tvector(k-1)+dT;
    
    
    if(k<=delayT)
        ToutputVector(k)=0;
    else
        ToutputVector(k)=Tvector(k-delayT);
    end
end

eNormalT = Tvector;
eNormalTout = ToutputVector;
eNormalV=hVector;

% figure
% plot(t,hVector)
% title("Napełnianie zbiornika" + newline + "symulacja metodą Eulera")
% xlabel("t[s]");
% ylabel("h[cm]")
% hold off
% 
% figure
% plot(t,Tvector)
% hold on
% plot(t,ToutputVector)
% title("Napełnianie zbiornika" + newline+"tempteratura"+newline + "symulacja metodą Eulera")
% xlabel("t[s]");
% ylabel("T[\circC]")
% legend("temperatura w zbiorniku", "temperatura na wyjsciu", 'Location','best')
% hold off


% Euler modified
for k=2:length(t)
    h=hVector(k-1);
    T=Tvector(k-1);

    if(k<delayC)
      delayFc=0;  
    else
       delayFc=1;
    end
    
   pomV=V(k-1)+0.5*Tp*dVdt(h,delayFc, Finputs);
   pomh=heightFromVolume(pomV);
   
   dV=Tp*dVdt(pomh,delayFc, Finputs);
   V(k)=V(k-1)+dV;
   hVector(k)=heightFromVolume(V(k));
   
   pomT=T + 0.5*Tp*dTdt(V(k),T,delayFc,Finputs,Tinputs);
   dT=Tp*dTdt(V(k-1),pomT,delayFc,Finputs,Tinputs);
   Tvector(k)=Tvector(k-1)+dT;
   
   if(k<=delayT)
       ToutputVector(k)=0;
   else
       ToutputVector(k)=Tvector(k-delayT);
   end
   
end

eModifiedT = Tvector;
eModifiedTout = ToutputVector;
eModifiedV = hVector;

% figure
% plot(t,hVector)
% title("Napełnianie zbiornika" + newline + "symulacja zmodyfikowaną metodą Eulera")
% xlabel("t[s]");
% ylabel("h[cm]")
% hold off

% figure
% plot(t,Tvector)
% title("Temperatura w zbiorniku"+newline+"podczas napełniania zbiornika")
% xlabel("t[s]");
% ylabel("T[\circC]")
% hold off

% rungy-kutta
for k=2:length(t)
    h=hVector(k-1);
    T=Tvector(k-1);

    if(k<delayC)
      delay=0;  
    else
       delay=1;
    end
    
    kV1= dVdt(h, delay, Finputs);
    kV2= dVdt(heightFromVolume(V(k-1) + Tp/2*kV1),delay,Finputs);
    kV3= dVdt(heightFromVolume(V(k-1) + Tp/2*kV2),delay,Finputs);
    kV4= dVdt(heightFromVolume(V(k-1) + Tp*kV3),delay,Finputs);
    
    
    dV=Tp/6*(kV1+2*kV2+2*kV3+kV4);
    V(k)=V(k-1)+dV;
    hVector(k)=heightFromVolume(V(k));
            
    kT1= dTdt(V(k-1),T,delay,Finputs,Tinputs);
    kT2= dTdt(V(k-1) + Tp*V(k-1)/2,T + Tp/2*kT1,delay,Finputs,Tinputs);
    kT3= dTdt(V(k-1) + Tp*V(k-1)/2,T + Tp/2*kT2,delay,Finputs,Tinputs);
    kT4= dTdt(V(k-1) + Tp*V(k-1),T + Tp*kT3,delay,Finputs,Tinputs);       
   
    dT=Tp/6*(kT1+2*kT2+2*kT3+kT4);
    Tvector(k)=Tvector(k-1)+ dT;

    if(k<=delayT)
        ToutputVector(k)=0;
    else
        ToutputVector(k)=Tvector(k-delayT);
    end
end

rkV = hVector;
rkT = Tvector;
rkTout = ToutputVector;

% figure
% plot(t,hVector)
% title("Napełnianie zbiornika" + newline + "symulacja metodą Rungego Kutty")
% xlabel("t[s]");
% ylabel("h[cm]")
% hold off

volumeFigure=figure
plot(t,eNormalV, 'r')
hold on
plot(t,eModifiedV,'g')
hold on
plot(t, rkV, 'b')
title("Napełnianie zbiornika" + newline + "symulacja metodą trzech różnych metod")
legend("zwykła metoda Eulera", "zmodyfikowana metoda Eulera", "metoda Rungego Kutty", 'Location', 'best');
xlabel("t[s]");
ylabel("h[cm]")
hold off

temperatureFigure=figure
plot(t,eNormalT, 'b')
hold on
plot(t,eNormalTout,'g')
hold on
plot(t,eModifiedT, 'r')
hold on
plot(t,eModifiedTout,'c')
hold on
plot(t,rkT, 'k')
hold on
plot(t,rkTout,'m')

title("Napełnianie zbiornika" + newline +"temperatura"+newline + "symulacja metodą trzech różnych metod")
legend("zwykła metoda Eulera T w zbiorniku","zwykła metoda Eulera T na wyjściu" ,...
    "zmodyfikowana metoda Eulera T w zbiorniku", "zmodyfikowana metoda Eulera T na wyjściu", ...
    "metoda Rungego Kutty T w zbiorniku","metoda Rungego Kutty T na wyjściu",...
    'Location', 'best');
xlabel("t[s]");
ylabel("T[\circC]")
hold off

fileName="object-simulation";
overLeafFilePath="img/object-simulation/";
path=overLeafFilePath;

name="method-comparison-volume";
caption="Symulacja napełniania zbironika przy użyciu trzech różnych metod obliczeń numerycznych";
label="fig:method-comparison-volume";

saveFigure(volumeFigure, path,name,fileName,overLeafFilePath,caption,label);


name="method-comparison-temperature";
caption="Symulacja zmian temperatury w zbiorniku i na wyjściu zbironika przy użyciu trzech różnych metod obliczeń numerycznych";
label="fig:method-comparison-volume";

saveFigure(temperatureFigure, path,name,fileName,overLeafFilePath,caption,label);
