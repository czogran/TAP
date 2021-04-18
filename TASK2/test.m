Kp = 0.1;
Ki = 0;   % No integrator
Kd = 0.5;
Tf = 0;
C = pid(Kp,Ki,Kd,Tf);

T=tf(1,[1 ,1]);


step(series(C,T));