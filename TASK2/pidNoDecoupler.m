transferModel;

% PID Controller Design at the Command Line

% s = tf('s','TimeUnit','minutes');
% G = [87.8 -86.4 ; 108.2 -109.6]/(75*s+1);
% G.InputName = {'L','V'};
% G.OutputName = 'y';

G=tf(transmit);

D = tunableGain('Decoupler',eye(2));

% D.u = {'e1','e2'};
D.u = 'e';

% D.y = {'pL','pV'};
D.y = {'V1','T1'};


C_Fh = tunablePID('C_V','pid');  
% C_L.TimeUnit = 'minutes';
C_Fh.u = 'e(1)'; 
C_Fh.y = 'Fh';

C_Fh.Kp.Value=1000;
C_Fh.Ki.Value=0;
C_Fh.Kd.Value=0;

C_Fc = tunablePID('C_T','pid');  
% C_V.TimeUnit = 'minutes';
C_Fc.u = 'e(2)'; 
C_Fc.y = 'Fc';


C_Fc.Kp.Value=-0.4;
C_Fc.Ki.Value=-0.003;


Sum1 = sumblk('e1 = r1 - V');
Sum2 = sumblk('e2 = r2 - T');
Sum = sumblk('e = r - y',2);



% CLry = connect(G,D,C_Fh,C_Fc,Sum1,'r1','V',Sum2,'r2','T');
% CLry = connect(G,D,C_Fh,C_Fc,Sum1,Sum2,{'r1','r2'},{'V', 'T'});
% CLry = connect(G,D,C_Fh,C_Fc,Sum,'r','y');
CLry = connect(G,C_Fh,C_Fc,Sum,'r','y',{'Fh','Fc'}, 'u');



[y t x]=lsim(CLry, ones(2,10000).*[81;33.57],1:10000);

figure
plot(y(:,1))
hold off

figure
plot(y(:,2))
hold off