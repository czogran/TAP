init

% consts
a=7;
C=0.7;

% sampe time
Tp=2; 

transmit = tf(1, [C,a* 0.5/sqrt(h0)],'InputDelay',180);

V0=volume(81);
transmit2 = tf(1, [1 ,a* 0.25/sqrt(sqrt(C*(V0^3)))],'InputDelay',180);
transmit3 = tf(1, [2*C*h0,a* 0.5/sqrt(h0) ],'InputDelay',180);


 Numerator = {[Th-T0] [0.8 72.3]}; %Numerators of u_1 and u_2
 Denominator = {[V0, (Fh+Fd+Fcin)] [1 21.8 60]}; %Denominators of u_1 and u_2
 transformT = tf(Numerator,Denominator)

figure
step(transformT)


method='zoh';
discreteTransmit=c2d(transmit,Tp,method); % transmitancja dyskretna

% figure
% step(transmit)
% hold on
% step(discreteTransmit)
% hold off

