%%--SYMULACJA DZIAŁANIA OBIEKTU--%
%Polecenie: Zasymulować działanie obiektu w Matlabie

%init data
init

t=0:1000;

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
plot(t,Tvector)
hold on
plot(t,ToutputVector)
title("Napełnianie zbiornika" + newline+"tempteratura"+newline + "symulacja metodą Eulera")
xlabel("t[s]");
ylabel("T[\circC]")
legend("temperatura w zbiorniku", "temperatura na wyjsciu", 'Location','best')
hold off