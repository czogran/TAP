% PID Controller Design at the Command Line

s = tf('s','TimeUnit','minutes');
G = [87.8 -86.4 ; 108.2 -109.6]/(75*s+1);
G.InputName = {'L','V'};
G.OutputName = 'y';

D = tunableGain('Decoupler',eye(2));
D.u = 'e';
D.y = {'pL','pV'};

C_L = tunablePID('C_L','pi');  
C_L.TimeUnit = 'minutes';
C_L.u = 'pL'; 
C_L.y = 'L';

C_L.Kp.Value=1;

C_V = tunablePID('C_V','pi');  
C_V.TimeUnit = 'minutes';
C_V.u = 'pV'; 
C_V.y = 'V';


C_V.Kp.Value=-10;


Sum = sumblk('e = r - y',2);

CLry = connect(G,D,C_L,C_V,Sum,'r','y');


[y t x]=lsim(CLry, ones(2,1000),1:1000);
