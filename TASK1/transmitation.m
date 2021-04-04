init

% consts
a=7;
C0=0.7;

% sampe time
Tp=1; 

addconst = 0;

%transmit = tf({1 1 1}, {[1, a/(4*sqrt(sqrt(h0^3*C0)))] [1, a/(4*sqrt(sqrt(h0^3*C0)))] [1, a/(4*sqrt(sqrt(h0^3*C0)))]}, 'InputDelay', [0; 120; 0]);
h0 = 81;
V0 = volume(h0);

if addconst == 1
    B = [1, 1, 1, 0, 0, 0; Th0/V0 - T0/V0, Tc0/V0 - T0/V0, Td0/V0 - T0/V0, Fh0/V0, Fc0/V0, Fd0/V0;0, 0, 0, 0, 0, 0];
    A = [-a*0.25/(sqrt(sqrt((V0^3)*C0))), 0, 1; Fh0*T0/(V0^2) + Fc0*T0/(V0^2) + Fd0*T0/(V0^2) - Fh0*Th0/(V0^2) - Fc0*Tc0/(V0^2) - Fd0*Td0/(V0^2), -(Fh0/V0 + Fc0/V0 + Fd0/V0), 0; 0, 0, 0];
    C = [1, 0, 0;0, 1, 0; 0, 0, 1];
    D = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0 ,0 , 0; 0, 0, 0, 0, 0, 0];
else
    B = [1, 1, 1, 0, 0, 0; Th0/V0 - T0/V0, Tc0/V0 - T0/V0, Td0/V0 - T0/V0, Fh0/V0, Fc0/V0, Fd0/V0];
    A = [-a*0.25/(sqrt(sqrt((V0^3)*C0))), 0; Fh0*T0/(V0^2) + Fc0*T0/(V0^2) + Fd0*T0/(V0^2) - Fh0*Th0/(V0^2) - Fc0*Tc0/(V0^2) - Fd0*Td0/(V0^2), -(Fh0/V0 + Fc0/V0 + Fd0/V0)];
    C = [1, 0;0, 1];
    D = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0 ,0 , 0];
end




transmit = ss(A, B, C, D, 'InputDelay', [0, 180, 0, 0, 0, 0]);


method='zoh';
%discreteTransmit=c2d(transmit,Tp,method); % transmitancja dyskretna
t = 1:2000;
u = [ones(1,2000)*Fh;ones(1,2000)*Fcin;ones(1,2000)*Fd;ones(1,2000)*Th;ones(1,2000)*Tc;ones(1,2000)*Td];

figure
if addconst == 1
    lsim(transmit, u, t, [V0, T0, -a/sqrt(sqrt(C0))*(sqrt(sqrt(V0))-0.25*V0/sqrt(sqrt(V0^3)))])
else
    lsim(transmit, u, t, [V0, T0])
end
hold on
%step(discreteTransmit)
hold off

