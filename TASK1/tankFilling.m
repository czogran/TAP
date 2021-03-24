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
Tinputs=[Fh,Fc,Fd];

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

figure
plot(t,hVector)
title("Napełnianie zbiornika" + newline + "symulacja metodą Eulera")
xlabel("t[s]");
ylabel("h[cm]")
hold off

figure
plot(t,Tvector)
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

eModified = hVector;

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
%     
%      k1= dVdt(h, delay, Finputs);
%     k2= dVdt(h+Tp/2,delay,Finputs);
%     k3= dVdt(h + Tp/2,delay,Finputs);
%     k4= dVdt(h+Tp,delay,Finputs);
% %     
%        k1= dVdt(h, delay, Finputs);
%     k2= dVdt(heightFromVolume(V(k-1) + Tp/2*k1),delay,Finputs);
%     k3= dVdt(heightFromVolume(V(k-1) + Tp/2*k2),delay,Finputs);
%     k4= dVdt(heightFromVolume(V(k-1) + Tp*k3),delay,Finputs);
    
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

name="method-comparison-flow";
fileName="object-simulation";
overLeafFilePath="img/object-simulation/";
caption="Symulacja napełniania zbironika przy użyciu trzech różnych metod obliczeń numerycznych";
label="fig:method-comparison-flow";

saveFigure(gcf, 'img/object-simulation/',name,fileName,overLeafFilePath,caption,label);
