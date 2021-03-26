init

% consts
a=7;
C=0.7;

% sampe time
Tp=2; 

transmit = tf(1, [C,a* 0.5/sqrt(h0)]);

method='zoh';
discreteTransmit=c2d(transmit,Tp,method); % transmitancja dyskretna

figure
step(transmit)
hold on
step(discreteTransmit)
hold off

