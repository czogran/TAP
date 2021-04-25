transferModel;

% PID Controller Design at the Command Line

G=tf(transmit);

PID_Fh = tunablePID('PID_Fh','pid');  
PID_Fh.u = 'e(1)'; 
PID_Fh.y = 'Fh';

PID_Fh.Kp.Value=1;
PID_Fh.Ki.Value=0;
PID_Fh.Kd.Value=0;

PID_Fc = tunablePID('PID_Fc','pid');  
PID_Fc.u = 'e(2)'; 
PID_Fc.y = 'Fc';


PID_Fc.Kp.Value=-0.03;
PID_Fc.Ki.Value=-0.005;

Sum = sumblk('e = r - y',2);

AP_Fh = AnalysisPoint('Fh');
% AP_Fc = AnalysisPoint('ch');
% 
% p=[PID_Fh 0;
%     0 PID_Fc;
%     0 0;
%     0 0;
%     0 0 ;
%     0 0];
% 
% T = feedback(G*AP_Fh*p,ones(2,2)');      % closed loop r->y


CLry = connect(G,PID_Fh,PID_Fc,Sum,'r','y',{'e','Fh','Fc'});

[y t x]=lsim(CLry, [ones(2,10000).*[V0;T0],ones(2,10000).*[V0-100;Th]],1:20000);


figure
plot(y(:,1))
hold off

figure
plot(y(:,2))
hold off


TuyE = getIOTransfer(CLry,'r(2)','e(2)');
Tuy = getIOTransfer(CLry,'r(1)','Fh');

e2=lsim(TuyE,y(:,2),1:length(y(:,1)));
U=lsim(Tuy,y(:,1),1:length(y(:,1)));
% U1=lsim(Tuy1,ones(1,10000)*81,1:10000);

figure
plot(U)
hold off