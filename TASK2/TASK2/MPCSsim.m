%%--SYMULACJA DZIAŁANIA OBIEKTU--%
%Polecenie: Zasymulować działanie obiektu w Matlabie

%init data and inports
step

N = size(S);
N = N(3);
Nu = N;
N = 5;
Nu = 4;

[M, MP] = DMCmatrices(S, N, Nu);
dimM = size(M);

psi = 1;
lambda = 1;

Psi = eye(N*2)*psi;
Lambda = eye(Nu*2)*lambda;

%Macierz K
K=(M'*Psi*M+Lambda)^(-1)*M';
dimK = size(K);

K11 = zeros(1, dimK(2)/2);
K12 = zeros(1, dimK(2)/2);
K21 = zeros(1, dimK(2)/2);
K22 = zeros(1, dimK(2)/2);

for i = 1:2:dimK(2)
    K11((i+1)/2) = K(1, i);
    K12((i+1)/2) = K(1, i+1);
    K21((i+1)/2) = K(2, i);
    K22((i+1)/2) = K(2, i+1);
end

Ke = [sum(K11), sum(K12); sum(K21), sum(K22)];

t = 1:100000;
%% tank filling
hVector=ones(length(t),1)*h0;

VVector = ones(length(t),1)*V0;
TVector = ones(length(t),1)*T0;
ToutputVector=TVector;

Yzad = ones(length(t), 2).*[84; 35];

Finputs=[Fh,Fc,Fd];
Tinputs=[Th,Tc,Td];
Tin = Tinputs;
Fin = Finputs;
delay = 1;
Tp = 1;
% rungy-kutta
for k=2:length(t)
    
    V = VVector(k-1);
    Tout = ToutputVector(k-1);
    T = TVector(k-1);
    
   
    
    dU = Ke * Yzad(k, :) + 
    
    
    
    
    
    
    
    if(k<delayC)
        Fin = [Fin(1), Fc, Fin(3)];
    end
    
    kV1= dVdt(heightFromVolume(V), delay, Fin);
    kV2= dVdt(heightFromVolume(V + Tp/2*kV1),delay,Fin);
    kV3= dVdt(heightFromVolume(V + Tp/2*kV2),delay,Fin);
    kV4= dVdt(heightFromVolume(V + Tp*kV3),delay,Fin);
    
    kT1= dTdt(V,T,delay,Finputs,Tinputs);
    kT2= dTdt(V + Tp/2*kV1,T + Tp/2*kT1,delay,Fin,Tin);
    kT3= dTdt(V + Tp/2*kV2,T + Tp/2*kT2,delay,Fin,Tin);
    kT4= dTdt(V + Tp*kV3,T + Tp*kT3,delay,Fin,Tin);
   
    dV=Tp/6*(kV1+2*kV2+2*kV3+kV4);
    VVector(k)=V+dV;

    
    dT=Tp/6*(kT1+2*kT2+2*kT3+kT4);   
    hVector(k)=heightFromVolume(VVector(k));
    %dT=dTdt(V(k),T,delayFc,Finputs,Tinputs);
    TVector(k)=TVector(k-1)+dT;
    if(k<=delayT)
        ToutputVector(k)=T0;
    else
        ToutputVector(k)=TVector(k-delayT);
    end
    
    %d = [VVector(k); ToutputVector(k)] - (A*[V;Tout] + B*[Fin, Tin]');
    d = [VVector(k); ToutputVector(k)] - [VVector(k-1); ToutputVector(k-1)] - (A*[V;Tout] + B*[Fin, Tin]');
%     Tvector(k)=Tvector(k-1)+ dVdTdt/V(k);
end

rkT = TVector;
rkTout = ToutputVector;

rk = hVector;


figure
plot(t, rk, 'b')
title("Napełnianie zbiornika" + newline + "symulacja metodą trzech różnych metod")
legend("zwykła metoda Eulera", "zmodyfikowana metoda Eulera", "metoda Rungego Kutty", 'Location', 'best');
xlabel("t[s]");
ylabel("h[cm]")
hold off

figure
plot(t,rkT)
title("Napełnianie zbiornika" + newline+"tempteratura"+newline + "symulacja metodą RK4")
xlabel("t[s]");
ylabel("T[\circC]")
legend("temperatura w zbiorniku", "temperatura zlinearyzowana", 'Location','best')
hold off
