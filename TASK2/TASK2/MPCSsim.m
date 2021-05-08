%%--SYMULACJA DZIAŁANIA OBIEKTU--%
%Polecenie: Zasymulować działanie obiektu w Matlabie

%init data and inports
step
sys = ss(A,B,C,D);
sysDys = c2d(sys,1);
A = sysDys.A;
B = sysDys.B;
B(1, :) = B(1,:)/ratio;
B = B(1:2, 1:2);

N = size(S);
N = N(3);
Nu = N;

N = 1500;
Nu = 1500;

N = 1500;
Nu = 1500;

M=zeros(N*2,Nu*2);
for i=1:N, M((i-1)*2+1:(i-1)*2+2,1:2)=S(:,:,min(i,Dyn)); end
for i=2:Nu
    M(:,(i-1)*2+1:(i-1)*2+2)=[zeros((i-1)*2,2); M(1:(N-i+1)*2,1:2)];
end
dimM = size(M);

psi = 1;
lambda = 0.5;

Psi = eye(N*2)*psi;
Lambda = eye(Nu*2)*lambda;

%Macierz K
K=((M'*Psi*M+Lambda)^(-1))*M'*Psi;
dimK = size(K);

K1 = K(1:2, :);

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

Vp = eye(size(A));

for i = 1:N-1
    acc = [0, 0;0, 0];
    for j = 1:i
        acc = A^j + acc;
    end
    Vp = [Vp; acc + eye(size(A))];
end
V0 = V0/ratio;

t = 1:1000;
%% tank filling
hVector=ones(length(t),1)*h0;

VVector = ones(length(t),1)*V0;
TVector = ones(length(t),1)*T0;
ToutputVector=TVector;

T0 = 33.57143;

Yzad = ones(length(t), 2).*[volume(70)/ratio, 40];


Finputs=[Fh,Fc,Fd];
Tinputs=[Th,Tc,Td];

FinVector = ones(length(t), 3).*Finputs;
TinVector = ones(length(t), 3).*Tinputs;

E = 0;

Tin = Tinputs;
Fin = Finputs;
delay = 1;
Tp = 1;
dist = [0; 0];
% rungy-kutta
for k=2:length(t)
    
    V = VVector(k-1);
    Tout = ToutputVector(k-1);
    T = TVector(k-1);
    
    acc = [0; 0];
    
    for i = 1:N
        acc = acc + K1(:, 2*i-1:2*i)*(C*((A^i)*[V - V0;Tout - T0]) + C*(Vp(2*i -1: 2*i, :))*(B*(FinVector(k-1, 1:2) - [Fh,Fc])') + dist);
    end
    
    
    dU = Ke * (Yzad(k, :) - [V0, T0])' - acc;
    
    FinVector(k, 1:2) = dU + FinVector(k-1, 1:2)';
    
    FinVector(k, 1) = min(FinVector(k, 1), 100); 
    FinVector(k, 2) = min(FinVector(k, 2), 100);
    
    FinVector(k, 1) = max(FinVector(k, 1), 0); 
    FinVector(k, 2) = max(FinVector(k, 2), 0);
    
    Fin = FinVector(k, :);
    Tin = TinVector(k, :);
    
    if(k<delayC)
        Fin = [Fin(1), Fc, Fin(3)];
    end
    
    kV1= dVdt(heightFromVolume(V*ratio), delay, Fin);
    kV2= dVdt(heightFromVolume(V*ratio + Tp/2*kV1),delay,Fin);
    kV3= dVdt(heightFromVolume(V*ratio + Tp/2*kV2),delay,Fin);
    kV4= dVdt(heightFromVolume(V*ratio + Tp*kV3),delay,Fin);
    
    kT1= dTdt(V*ratio,T,delay,Fin,Tin);
    kT2= dTdt(V*ratio + Tp/2*kV1,T + Tp/2*kT1,delay,Fin,Tin);
    kT3= dTdt(V*ratio + Tp/2*kV2,T + Tp/2*kT2,delay,Fin,Tin);
    kT4= dTdt(V*ratio + Tp*kV3,T + Tp*kT3,delay,Fin,Tin);
   
    dV=Tp/6*(kV1+2*kV2+2*kV3+kV4)/ratio;
    VVector(k)=V+dV;

    
    dT=Tp/6*(kT1+2*kT2+2*kT3+kT4);   
    hVector(k)=heightFromVolume(VVector(k)*ratio);
    %dT=dTdt(V(k),T,delayFc,Finputs,Tinputs);
    TVector(k)=TVector(k-1)+dT;
    if(k<=delayT)
        ToutputVector(k)=T0;
    else
        ToutputVector(k)=TVector(k-delayT);
    end
    
    %d = [VVector(k); ToutputVector(k)] - (A*[V;Tout] + B*[Fin, Tin]');
    dist = [VVector(k); ToutputVector(max(k, 1))] - (A*[V;ToutputVector(max(k - 1, 1))] + B*FinVector(max(k - 1, 1), 1:2)');
    
    E = E + sum((Yzad(k, :) - [VVector(k), TVector(k)]).^2);
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
plot(t,rkTout)
title("Napełnianie zbiornika" + newline+"tempteratura"+newline + "symulacja metodą RK4")
xlabel("t[s]");
ylabel("T[\circC]")
legend("temperatura w zbiorniku", "temperatura zlinearyzowana", 'Location','best')
hold off

figure
plot(FinVector(:,1))

figure
plot(FinVector(:,2))
