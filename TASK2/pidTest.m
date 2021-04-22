%% RÓŻNICA TAKA ŻE ZAMIAST PID TUNNABLE ZWYKŁY PID, może jest jakas roznica...

transferModel

H = pid(1,0.01); 
H.InputName = 'e(1)';  
H.OutputName = 'Fh';

C = pid(-1,-0.005); 
C.InputName = 'e(2)';  
C.OutputName = 'Fc';

G=tf(transmit);

Sum = sumblk('e = r - y',2);

CLry = connect(G,H,C,Sum,'r','y',{'Fh','Fc'});



[y t x]=lsim(CLry, ones(2,10000).*[81;33.57],1:10000);

figure
plot(y(:,1))
hold off

figure
plot(y(:,2))
hold off
