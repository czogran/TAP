transferModelNoDelays;

UseBackCalculation = true; % enable/disable anti-windup (back calculation)
G=tf(transmit);
KpFh = 0.1; KiFh = 0.001; % PI controller gains (parallel)
KpFc = -20; KiFc = -0.5; % PI controller gains (parallel)

Ts = 1; % controller sample time
tau = 1; % reset time constant
UB = 300; LB = 0; % saturation limits

actionIFh = 0;
actionIFc = 0;

% closed-loop simulation (200 steps)
N = 4000;

% initial condition
r1 = [h0, T0];
r2=[h0-10, T0+5];
rVector=[ones(N-2000,2).*r1;[ones(2000,2).*r2]];

% vector pre allocation
y = [0,0];

offset=max(delayC,delayT);
yVec = zeros(N+offset,2);
uVecFh= zeros(N+offset,1);
uVecFc= zeros(N+offset,2);

for ct=1:N
    % error
    e = rVector(ct,:) - y;
    
    % control action
    actionPFh = KpFh*e(1);
    actionPFc = KpFc*e(2);

    uFh = [actionPFh + actionIFh];
    uFc = [actionPFc + actionIFc];

    % saturation control action
    uFhSat = max(min(uFh,UB),LB);
    uFcSat = max(min(uFc,UB),LB);

    u=[uFhSat, uFcSat];
    
    % anti windup
    if UseBackCalculation
        actionIFh = actionIFh + (KiFh*e(1) + (uFhSat-uFh)/tau)*Ts;
        actionIFc = actionIFc + (KiFc*e(2) + (uFcSat-uFc)/tau)*Ts;
    else
        actionIFh = actionIFh + KiFh*e(1)*Ts;
        actionIFc = actionIFc + KiFc*e(2)*Ts;
    end
    
    
    V=volume(y(1));
    Finputs=[u,Fd];
    Tinputs=[Th,Tc,Td];
    
    if(ct*Ts<delayC)
      delay=0;  
    else
       delay=1;
    end
    
    kV1= dVdt(heightFromVolume(V), delay, Finputs);
    kV2= dVdt(heightFromVolume(V + Ts/2*kV1),delay,Finputs);
    kV3= dVdt(heightFromVolume(V + Ts/2*kV2),delay,Finputs);
    kV4= dVdt(heightFromVolume(V + Ts*kV3),delay,Finputs);
    
    dV=Ts/6*(kV1+2*kV2+2*kV3+kV4);
    V=V+dV;
    y(1)=heightFromVolume(V);
            
    T=y(2);
    kT1= dTdt(V,T,delay,Finputs,Tinputs);
    kT2= dTdt(V + Ts*V/2,T + Ts/2*kT1,delay,Finputs,Tinputs);
    kT3= dTdt(V + Ts*V/2,T + Ts/2*kT2,delay,Finputs,Tinputs);
    kT4= dTdt(V + Ts*V,T + Ts*kT3,delay,Finputs,Tinputs);       
    
    dT=Ts/6*(kT1+2*kT2+2*kT3+kT4);
    y(2)=T+ dT;
    
    yVec(ct,1) = y(1);
    yVec(ct+delayT,2) = y(2);
    uVecFh(ct) = u(1,1);
    uVecFc(ct+delayC) = u(1,2);
end

% remove offset
yVec = yVec(1:N,:);
uVecFh= uVecFh(1:N);
uVecFc= uVecFc(1:N);

clf;
figure(1);
plot(1:N,yVec(:,1),'.r');
hold on
plot(rVector(:,1),'--b');
xlabel('Time'); ylabel('Signal');
if UseBackCalculation
    title('With Anti-Windup (Back Calculation)');
else
    title('Without Anti-Windup');
end
hold off

figure(2);
plot(1:N,yVec(:,2),'.r');
hold on
plot(rVector(:,2),'--b');

% plot(1:N,y_vec(:,2),'*r',1:N,u_vec,'+b');
xlabel('Time'); ylabel('Signal');
legend('Model Response','Control Signal')
if UseBackCalculation
    title('With Anti-Windup (Back Calculation)');
else
    title('Without Anti-Windup');
end


figure
plot(uVecFh,'--r')
hold on
plot(uVecFc,'--g')
title("u")
legend("Fh","Fc")
hold off






% 
% transferModelNoDelays;
% 
% G=tf(transmit);
% 
% PID_Fh = tunablePID('PID_Fh','pid');  
% PID_Fh.u = 'e(1)'; 
% PID_Fh.y = 'Fh';
% 
% PID_Fh.Kp.Value=2;
% PID_Fh.Ki.Value=10;
% PID_Fh.Kd.Value=0.001;
% 
% PID_Fc = tunablePID('PID_Fc','pid');  
% PID_Fc.u = 'e(2)'; 
% PID_Fc.y = 'Fc';
% 
% PID_Fc.Kp.Value=0;
% PID_Fc.Ki.Value=-150;
% PID_Fc.Kd.Value=0.001;
% 
% Sum = sumblk('e = r - y',2);
% 
% CLry = connect(PID_Fh,PID_Fc,Sum,{'r','y'},{'Fh', 'Fc'});
% 
% % CLry = connect(G,PID_Fh,PID_Fc,Sum,'r','y',{'e','Fh','Fc'});
% tSim=500;
% offset= 200;
% 
% r=[ones(round(tSim/10),2).*[V0-1000, T0];ones(round(tSim/1),2).*[V0, T0+20];ones(round(tSim/1),2).*[V0, T0];ones(tSim,2).*[V0, T0+10];];
% 
% r=ones(tSim+offset,2).*[V0, T0];
% 
% Fh = ones(offset+tSim, 1)*Fh0;
% Fc = ones(offset+tSim, 1)*Fc0;
% 
% y=ones(tSim+offset,2).*[V0, T0];
% 
% % ,'Fd','Th','Tc','Td'
% constControls=ones(101,1)*[Fd0,Th0,Tc0,Td0];
% 
% c={};
% x=[0,0 ,0 ,0];
% x2=[0,0,0];
% 
% for k=1:tSim
% %         controls = lsim(PID_Fh,[r(k:k+1,1)-y(k:k+1,1)],0:1);
% 
%     [controls, t, x] = lsim(CLry,[r(k:k+1,:),y(k:k+1,:)],0:1,x);
%     x=x(end,:)';
%     c{k}=controls;
%     
%     if controls(2,1)<0
%         controls(2,1)=0;
%     end
%     
%     if controls(2,2)<0
%         controls(2,2)=0;
%     end
% %     
%     Fh(k+1)=controls(2,1);
%     Fc(k+1+delayC) = controls(2,2);
% %     
%     [ySim t1 x2] =lsim(G,[Fh(k:k+100),Fc(k:k+100),constControls],0:100);
% %     x2 = x2(end,:)';
%     y(k+1,1)=ySim(2,1);
% 
%     y(k+delayT,2)=ySim(1,2);
% end
% % [y t x]=lsim(CLry, [ones(2,10000).*[V0;T0],ones(2,10000).*[V0-100;Th]],1:20000);
% 
% figure
% plot(y(2:tSim,1))
% title("V")
% hold off
% 
% figure
% plot(y(2:tSim,2))
% title("T")
% hold off
% 
% figure
% plot(Fh)
% title("Fh")
% hold off
% 
% figure
% plot(Fc)
% title("Fc")
% hold off