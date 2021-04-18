%% DECOUPLING

% higher power first -> 10s^2+5s+3-> [10 5 3]
G11=exp(-5*s)/(24*s^2+10*s+4);
G12=exp(-7*s)*0.4*s/((1+9*s));

w=-G12/G11

G11= tf([1],[24 10 1], 'InputDelay',5);
G12 =tf([0.4],[9 1], 'InputDelay',7);
G21=tf([0.5],[5 1], 'InputDelay',6);
G22=tf([-0.4 2],[4 1], 'InputDelay',7);

% G = [G11, G12;
%     G21, G22];


% D12=-G12/G11;
% D22=-G21/G22;

K=[ 1 0.4;
    0.5 2];

K.*(K^-1)'
