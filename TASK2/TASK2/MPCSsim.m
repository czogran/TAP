%%--SYMULACJA DZIAŁANIA OBIEKTU--%
step
Umin = [0, 0];
Umax = [100, 100];
dUmin = [-0.2, -0.2];
dUmax = [0.2, 0.2];

sys = ss(A,B,C,D);
sysDys = c2d(sys,1);
A = sysDys.A;
B = sysDys.B;
B(1, :) = B(1,:)/ratio;
B = B(1:2, 1:2);

N = 500;
Nu = 60;

M=zeros(N*2,Nu*2);
for i=1:N, M((i-1)*2+1:(i-1)*2+2,1:2)=S(:,:,min(i,Dyn)); end
for i=2:Nu
    M(:,(i-1)*2+1:(i-1)*2+2)=[zeros((i-1)*2,2); M(1:(N-i+1)*2,1:2)];
end
dimM = size(M);

psi = 1;
lambda = 0.4;

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

t = 1:20000;
%% tank filling
hVector=ones(length(t),1)*h0;

VVector = ones(length(t),1)*V0;
TVector = ones(length(t),1)*T0;
ToutputVector=TVector;

T0 = 33.57143;

global Yzad

Yzad = ones(length(t), 2);
Yzad(1:2000, :) = Yzad(1:2000, :).*[volume(81)/ratio, 33.57];
Yzad(2001:5000, :) = Yzad(2001:5000, :).*[volume(90)/ratio, 38];
Yzad(5001:10000, :) = Yzad(5001:10000, :).*[volume(70)/ratio, 32];
Yzad(10001:15000, :) = Yzad(10001:15000, :).*[volume(70)/ratio, 40];
Yzad(15001:20000, :) = Yzad(15001:20000, :).*[volume(60)/ratio, 33];

Finputs=[Fh,Fc,Fd];
Tinputs=[Th,Tc,Td];

FinVector = ones(length(t), 3).*Finputs;
TinVector = ones(length(t), 3).*Tinputs;

FinVector(1:4000, 3) = Fd;
FinVector(4001:8000, 3) = Fd + 5;
FinVector(8001:12000, 3) = Fd;
FinVector(12001:16000, 3) = Fd - 5;
FinVector(16001:20000, 3) = Fd + 5;
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
    
    dU(1) = max(dU(1), dUmin(1));
    dU(2) = max(dU(2), dUmin(1));
    
    dU(1) = min(dU(1), dUmax(1));
    dU(2) = min(dU(2), dUmax(1));
    
    FinVector(k, 1:2) = dU + FinVector(k-1, 1:2)';
    
    FinVector(k, 1) = min(FinVector(k, 1), Umax(1)); 
    FinVector(k, 2) = min(FinVector(k, 2), Umax(2));
    
    FinVector(k, 1) = max(FinVector(k, 1), Umin(1)); 
    FinVector(k, 2) = max(FinVector(k, 2), Umin(2));
    
    Fin = FinVector(k, :);
    Tin = TinVector(k, :);
    
    if(k<=delayC)
        Fin = [Fin(1), Fc, Fin(3)];
    else
        Fin = [Fin(1), FinVector(k - delayC, 2), Fin(3)];
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
    dist = [VVector(k) - V0; ToutputVector(max(k, 1)) - T0] - (A*[V - V0;ToutputVector(max(k - 1, 1)) - T0] + B*(FinVector(max(k - 1, 1), 1:2)' - [Fh; Fc]));
    
    E = E + sum((Yzad(k, :) - [VVector(k), TVector(k)]).^2);
%     Tvector(k)=Tvector(k-1)+ dVdTdt/V(k);
end
lorg = lambda;
lambda = floor(100*lambda);

u1 = FinVector(:,1);
u2 = FinVector(:,2);
Yzad1 = heightFromVolume(Yzad(:, 1).*ratio);
Yzad2 = Yzad(:, 2);

fileName = "distrkN" + string(N) + "Nu" + string(Nu) + "l" + string(lambda);
overLeafFilePath = "img/MPCSanaRK/";

heightFigure=figure;
plot(hVector,'r');
hold on
plot(Yzad1,'--b');
xlabel('Time'); ylabel('Signal');
title("Regulator MPCS"+newline+"h[cm]" );
xlabel("t[s]")
ylabel("h[cm]")
legend("h[cm]", "trajektoria zadana", 'Location','best')
hold off

tempFigure=figure;
plot(TVector,'.r');
hold on
plot(Yzad2,'--b');
title("Regulator MPCS"+newline+"T[\circC]");
xlabel('t[s]'); ylabel('T[\circC]');
legend('temperatura [\circC]','trajektoria zadana', 'Location','best')
hold off

controlPIFigure=figure;
plot(u1,'r')
hold on
plot(u2,'g')
plot(FinVector(:,3),'b')
title("u")
title("Regulator MPCS"+newline+"sterowanie u");
legend("Fh","Fcin","Fd", 'Location','best')
xlabel("t[s]")
ylabel("u[$\frac{cm^3}{s}$]",'Interpreter','latex')
hold off
caption="Wykresy dla regulatora MPCS, obiekt nieliniowy.";
label="fig:distMPCSRK"+'N' + string(N) + 'Nu' + string(Nu) + 'l' + string(lambda);

heightName="distMPCSRKH"+"N" + string(N) + "Nu" + string(Nu) + "l" + string(lambda);
tempName="distMPCSRKT"+"N" + string(N) + "Nu" + string(Nu) + "l" + string(lambda);
controlPIName="distMPCSRKControl"+"N" + string(N) + "Nu" + string(Nu) + "l" + string(lambda);

saveFiguresInColumn([heightFigure,tempFigure,controlPIFigure], "img/MPCSanaRK/",[heightName,tempName,controlPIName],fileName,overLeafFilePath,caption,label);


%figure
%plot(t, hVector, 'b')
%hold on
%plot(t, Yzad1, 'LineStyle', '--')
%title("Napełnianie zbiornika")
%xlabel("t[s]");
%ylabel("h[cm]")
%hold off

%figure
%plot(t,ToutputVector)
%hold on
%plot(t, Yzad2, 'LineStyle', '--')
%title("Tempteratura w zbiorniku")
%xlabel("t[s]");
%ylabel("T[\circC]")
%hold off

%figure
%plot(FinVector(:,1))

%figure
%plot(FinVector(:,2))