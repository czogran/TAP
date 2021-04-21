init

% consts
a=7;
C=0.7;

% sampe time
Tp=2; 

V0=volume(81);

transformVFc = tf(1, [1 ,a* 0.25/sqrt(sqrt(C*(V0^3)))],'InputDelay',180);
transformHFc = tf(1, [2*C*h0,a* 0.5/sqrt(h0) ],'InputDelay',delayC);

transformVFh = tf(1, [1 ,a* 0.25/sqrt(sqrt(C*(V0^3)))]);
transformHFh = tf(1, [2*C*h0,a* 0.5/sqrt(h0) ]);

transformVFd = tf(1, [1 ,a* 0.25/sqrt(sqrt(C*(V0^3)))]);
transformHFd = tf(1, [2*C*h0,a* 0.5/sqrt(h0) ]);




% TEMPERATURE
transformTFh=tf([Th-T0],[V0, (Fh+Fd+Fcin)], 'OutputDelay', delayT );
transformTFc=tf([Tc-T0],[V0, (Fh+Fd+Fcin)],'InputDelay',delayC, 'OutputDelay', delayT);
transformTFd=tf([Td-T0],[V0, (Fh+Fd+Fcin)], 'OutputDelay', delayT );

transformTV=tf ([((Fh+Fcin+Fd)*T0-(Fh*Th+Fcin*Tc+Fd*Td))/V0],[V0, (Fh+Fd+Fcin)], 'OutputDelay', delayT );

transformTTh=tf(Fh,[V0, (Fh+Fd+Fcin)], 'OutputDelay', delayT );
transformTTc=tf(Fcin,[V0, (Fh+Fd+Fcin)], 'InputDelay',delayC, 'OutputDelay', delayT );
transformTTd=tf(Fd,[V0, (Fh+Fd+Fcin)], 'OutputDelay', delayT );

% plotTransforms




% figure
% step(transmit)
% hold on
% step(discreteTransmit)
% hold off

