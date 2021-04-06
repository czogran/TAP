init

% consts
a=7;
C=0.7;

% sampe time
Tp=2; 

V0=volume(81);

transmitVFc = tf(1, [1 ,a* 0.25/sqrt(sqrt(C*(V0^3)))],'InputDelay',180);
transmitHFc = tf(1, [2*C*h0,a* 0.5/sqrt(h0) ],'InputDelay',180);

transmitVFh = tf(1, [1 ,a* 0.25/sqrt(sqrt(C*(V0^3)))]);
transmitHFh = tf(1, [2*C*h0,a* 0.5/sqrt(h0) ]);

transmitVFd = tf(1, [1 ,a* 0.25/sqrt(sqrt(C*(V0^3)))]);
transmitHFd = tf(1, [2*C*h0,a* 0.5/sqrt(h0) ]);

transmitVFcFigure=figure
step(transmitVFc)
title("Transmitancja objętości w zbiorniku po skoku Fc")

transmitVFhFigure=figure
step(transmitVFh)
title("Transmitancja objętości w zbiorniku po skoku Fh")

transmitVFdFigure=figure
step(transmitVFd)
title("Transmitancja objętości w zbiorniku po skoku Fd")







transmitHFcFigure=figure
step(transmitHFc)
title("Transmitancja wysokości słupa cieczy w zbiorniku po skoku Fc")

transmitHFhFigure=figure
step(transmitHFh)
title("Transmitancja wysokości słupa cieczy w zbiorniku po skoku Fh")

transmitHFhFigure=figure
step(transmitHFd)
title("Transmitancja wysokości słupa cieczy w zbiorniku po skoku Fd")


% figure
% step(transmitH)


 Numerator = {[Th-T0] [0.8 72.3]}; %Numerators of u_1 and u_2
 Denominator = {[V0, (Fh+Fd+Fcin)] [1 21.8 60]}; %Denominators of u_1 and u_2
 transformT = tf(Numerator,Denominator)

figure
stepplot(transformT)
title("ww")

% method='zoh';
% discreteTransmit=c2d(transmit,Tp,method); % transmitancja dyskretna
% 


% figure
% step(transmit)
% hold on
% step(discreteTransmit)
% hold off

