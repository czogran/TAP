init
a=7;
C0=0.7;

Dyn = 1500;
S = zeros(2, 2, Dyn);

% sampe time
Tp=1; 

addconst = 0;
Vconst = -a/sqrt(sqrt(C0))*(sqrt(sqrt(V0))-0.25*V0/sqrt(sqrt(V0^3)));
%transmit = tf({1 1 1}, {[1, a/(4*sqrt(sqrt(h0^3*C0)))] [1, a/(4*sqrt(sqrt(h0^3*C0)))] [1, a/(4*sqrt(sqrt(h0^3*C0)))]}, 'InputDelay', [0; 120; 0]);


if addconst == 1
    B = [1, 1, 1, 0, 0, 0; Th0/V0 - T0/V0, Tc0/V0 - T0/V0, Td0/V0 - T0/V0, Fh0/V0, Fc0/V0, Fd0/V0;0, 0, 0, 0, 0, 0];
    A = [-a*0.25/(sqrt(sqrt((V0^3)*C0))), 0, 1; Fh0*T0/(V0^2) + Fc0*T0/(V0^2) + Fd0*T0/(V0^2) - Fh0*Th0/(V0^2) - Fc0*Tc0/(V0^2) - Fd0*Td0/(V0^2), -(Fh0/V0 + Fc0/V0 + Fd0/V0), 0; 0, 0, 0];
    C = [1, 0, 0;0, 1, 0];
    D = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0 ,0 , 0];
else
    B = [1, 1, 1, 0, 0, 0; Th0/V0 - T0/V0, Tc0/V0 - T0/V0, Td0/V0 - T0/V0, Fh0/V0, Fc0/V0, Fd0/V0];
    A = [-a*0.25/(sqrt(sqrt((V0^3)*C0))), 0; Fh0*T0/(V0^2) + Fc0*T0/(V0^2) + Fd0*T0/(V0^2) - Fh0*Th0/(V0^2) - Fc0*Tc0/(V0^2) - Fd0*Td0/(V0^2), -(Fh0/V0 + Fc0/V0 + Fd0/V0)];
    C = [1, 0;0, 1];
    D = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0 ,0 , 0];
end

hVector=zeros(length(t), 2);

VVector=[];
VVector = zeros(length(t), 2);
TVector= zeros(length(t), 2);
ToutputVector=zeros(length(t), 2);

Fc=Fcin;
Finputs=[Fh,Fc,Fd];
Tinputs=[Th,Tc,Td];

stepStart = 5000;

FinVector(:, :, 1) = ones(length(t), 3).*Finputs;
FinVector(:, :, 2) = ones(length(t), 3).*Finputs;
FinVector(1, :, 1) = [0, 0, 0];
FinVector(1, :, 2) = [0, 0, 0];
FinVector(stepStart:length(t) - 1, :, 1) = ones((length(t) + 1)/2,3).*[FinVector(1000, 1) + 1, FinVector(1000, 2), FinVector(1000, 3)];
FinVector(stepStart:length(t) - 1, :, 2) = ones((length(t) + 1)/2,3).*[FinVector(1000, 1), FinVector(1000, 2) + 1, FinVector(1000, 3)];

TinVector = ones(length(t), 3).*Tinputs;
TinVector(1, :) =  [0, 0, 0];



% figure
% plot(t,Tvector)
% title("Temperatura w zbiorniku"+newline+"podczas napełniania zbiornika")
% xlabel("t[s]");
% ylabel("T[\circC]")
% hold off

doVol = 1;
hPrev = 0;

%% sim
for k=2:length(t)
    h1=hVector(k-1, 1);
    T1=TVector(k-1, 1);
    V1=VVector(k-1, 1);
    h2=hVector(k-1, 2);
    T2=TVector(k-1, 2);
    V2=VVector(k-1, 2);
   
    
    Fin1 = [FinVector(k, 1, 1), FinVector(max(k - delayC, 1), 2, 1), FinVector(k, 3, 1)];
    Fin2 = [FinVector(k, 1, 2), FinVector(max(k - delayC, 1), 2, 2), FinVector(k, 3, 2)];
    Tin = [TinVector(k, 1), TinVector(max(k - delayC, 1), 2), TinVector(k, 3)];
    VVector(k, 1)=V1 + A(1,:) * [V1, T1]' + B(1,:) * [Fin1, Tin]' + Vconst;
    VVector(k, 2)=V2 + A(1,:) * [V2, T2]' + B(1,:) * [Fin2, Tin]' + Vconst;
    
    
    hVector(k, 1)=heightFromVolume(VVector(k, 1));
    hVector(k, 2)=heightFromVolume(VVector(k, 2));
    %dT=dTdt(V(k),T,delayFc,Finputs,Tinputs);
    TVector(k, 1) = T1 + A(2,:) * [V1, T1]' + B(2,:) * [Fin1, Tin]';
    TVector(k, 2) = T2 + A(2,:) * [V2, T2]' + B(2,:) * [Fin2, Tin]';
    
    if(k<=delayT)
        ToutputVector(k, 1)=0;
        ToutputVector(k, 2)=0;
    else
        ToutputVector(k, 1)=TVector(k-delayT, 1);
        ToutputVector(k, 2)=TVector(k-delayT, 2);
    end
%     Tvector(k)=Tvector(k-1)+ dVdTdt/V(k);
end

S(1, 1, :) = (VVector(stepStart:stepStart+Dyn-1, 1) - V0)/ratio;
S(2, 1, :) = ToutputVector(stepStart:stepStart+Dyn-1, 1) - T0;
S(1, 2, :) = (VVector(stepStart:stepStart+Dyn-1, 2) - V0)/ratio;
S(2, 2, :) = ToutputVector(stepStart:stepStart+Dyn-1, 2) - T0;
s11 = S(1, 1, :);
s11 = s11(:);
s21 = S(2, 1, :);
s21 = s21(:);
s12 = S(1, 2, :);
s12 = s12(:);
s22 = S(2, 2, :);
s22 = s22(:);
rkT = TVector;
rkTout = ToutputVector;

rk = hVector;
