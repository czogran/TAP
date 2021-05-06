function E = MPCSEval(N, Nu, lambda, S, Yzad, FdVector)
%%--SYMULACJA DZIAŁANIA OBIEKTU--%
%Polecenie: Zasymulować działanie obiektu w Matlabie

%init data and inports
Th=62;
Fh=14;
%Cold water
Tc=23;
Fc=37;
%Discruption
Td=33;
Fd=12;

Th0=62;
Fh0=14;
%Cold water
Tc0=23;
Fc0=37;
%Discruption
Td0=33;
Fd0=12;

a = 7;
C0 = 0.7;
Dyn = 1500;
%OUTPUT
Tout=0;
F=0;

h=81;
h0=81;

V0 = volume(h0);


T=33.57143;
T0 = T;
%ADJUSTABLE SIZES
% h;
% Tout;

%CONTROL VALUES
% Fh;
% Fcin-> Fc delayed;

%VARIABLES
% sample time
Tp=1;
t=0:Tp:10000;

%DELAY
delayT=120;
delayC=180;
ratio = 80;

B = [1, 1, 1, 0, 0, 0; Th0/V0 - T0/V0, Tc0/V0 - T0/V0, Td0/V0 - T0/V0, Fh0/V0, Fc0/V0, Fd0/V0];
A = [-a*0.25/(sqrt(sqrt((V0^3)*C0))), 0; Fh0*T0/(V0^2) + Fc0*T0/(V0^2) + Fd0*T0/(V0^2) - Fh0*Th0/(V0^2) - Fc0*Tc0/(V0^2) - Fd0*Td0/(V0^2), -(Fh0/V0 + Fc0/V0 + Fd0/V0)];
C = [1, 0;0, 1];
A = A + eye(size(A));
B(1, :) = B(1,:)/ratio;
B = B(1:2, 1:2);

M=zeros(N*2,Nu*2);
for i=1:N, M((i-1)*2+1:(i-1)*2+2,1:2)=S(:,:,min(i,Dyn)); end
for i=2:Nu
    M(:,(i-1)*2+1:(i-1)*2+2)=[zeros((i-1)*2,2); M(1:(N-i+1)*2,1:2)];
end

psi = 1;

Psi = eye(N*2)*psi;
Lambda = eye(Nu*2)*lambda;

%Macierz K
K=((M'*Psi*M+Lambda)^(-1))*M';
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

temp = size(Yzad);

t = 1:temp(1);
%% tank filling
hVector=ones(length(t),1)*h0;

VVector = ones(length(t),1)*V0;
TVector = ones(length(t),1)*T0;
ToutputVector=TVector;

T0 = 33.57143;

Finputs=[Fh,Fc,Fd];
Tinputs=[Th,Tc,Td];

FinVector = ones(length(t), 3).*Finputs;
TinVector = ones(length(t), 3).*Tinputs;

if length(FdVector) == 1
    FinVector(:, 3) = ones(temp(1), 1) * FdVector;
else
    FinVector(:, 3) = FdVector;
end

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
    dist = [VVector(k); ToutputVector(max(k, 1))] - (A*[V;ToutputVector(max(k, 1))] + B*[Fin(1), FinVector(max(k - 1, 1), 2)]');
    
    E = E + sum((Yzad(k, :) - [VVector(k), TVector(k)]).^2);
%     Tvector(k)=Tvector(k-1)+ dVdTdt/V(k);
end

end
