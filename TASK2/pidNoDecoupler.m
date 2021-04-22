transferModel;

% PID Controller Design at the Command Line

G=tf(transmit);

PID_Fh = tunablePID('C_V','pid');  
PID_Fh.u = 'e(1)'; 
PID_Fh.y = 'Fh';

PID_Fh.Kp.Value=1000;
PID_Fh.Ki.Value=0;
PID_Fh.Kd.Value=0;

PID_Fc = tunablePID('C_T','pid');  
PID_Fc.u = 'e(2)'; 
PID_Fc.y = 'Fc';


PID_Fc.Kp.Value=-0.4;
PID_Fc.Ki.Value=-0.003;

Sum = sumblk('e = r - y',2);

CLry = connect(G,PID_Fh,PID_Fc,Sum,'r','y');

[y t x]=lsim(CLry, ones(2,10000).*[81;23.57],1:10000);

figure
plot(y(:,1))
hold off

figure
plot(y(:,2))
hold off