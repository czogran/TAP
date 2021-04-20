% PID Controller Design at the Command Line

sys = zpk([],[-1 -1 -1],1); 
sys = tf(1, [1 3 3 1]);
sys=transmit
[C_pi,info] = pidtune(sys,'PI')

T_pi = feedback(C_pi*sys, 1);
step(T_pi)

[C_pi_fast,info] = pidtune(sys,'PI',1.0)

T_pi_fast = feedback(C_pi_fast*sys,1);
step(T_pi,T_pi_fast)
axis([0 30 0 1.4])
legend('PI','PI,fast')

[C_pidf_fast,info] = pidtune(sys,'PIDF',1.0)

T_pidf_fast =  feedback(C_pidf_fast*sys,1);
step(T_pi_fast, T_pidf_fast);
axis([0 30 0 1.4]);
legend('PI,fast','PIDF,fast');

S_pi_fast = feedback(sys,C_pi_fast);
S_pidf_fast = feedback(sys,C_pidf_fast);
step(S_pi_fast,S_pidf_fast);
axis([0 50 0 0.4]);
legend('PI,fast','PIDF,fast');
