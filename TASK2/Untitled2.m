transferModelNoDelays;

UseBackCalculation = true; % enable/disable anti-windup (back calculation)
G=tf(transmit);
KpFh = 0.1; KiFh = 0.001; % PI controller gains (parallel)
KpFc = -20; KiFc = -0.5; % PI controller gains (parallel)

PID_Fh = tunablePID('PID_Fh','pi');  
PID_Fh.u = 'e(1)'; PID_Fh.y = 'Fh';
PID_Fh.Kp.Value=KpFh; PID_Fh.Ki.Value=0;

PID_FhI = tunablePID('PID_FhI','pi');  
PID_FhI.u = 'e(1)'; PID_FhI.y = 'Fh';
PID_FhI.Kp.Value=0; PID_FhI.Ki.Value=KiFh;


PID_Fc = tunablePID('PID_Fc','pi');  
PID_Fc.u = 'e(2)'; PID_Fc.y = 'Fc';
PID_Fc.Kp.Value=KpFc; PID_Fc.Ki.Value=0;

PID_FcI = tunablePID('PID_FcI','pi');  
PID_FcI.u = 'e(2)'; PID_FcI.y = 'Fc';
PID_FcI.Kp.Value=0; PID_FcI.Ki.Value=KiFc;


Sum = sumblk('e = r - y',2);

CLry = connect(PID_Fh,PID_Fc,Sum,{'r','y'},{'Fh', 'Fc'});
CLryI = connect(PID_FhI,PID_FcI,Sum,{'r','y'},{'Fh', 'Fc'});


Ts = 1; % controller sample time
tau = 1; % reset time constant
UB = 300; LB = 0; % saturation limits

actionIFh = 0;
actionIFc = 0;

actionIFh2 = 0;
actionIFc2 = 0;

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

constControls=ones(2,1)*[Fd0,Th0,Tc0,Td0];

x=[];
xI=[0 0];

cVector=[];

% D12=tf([0.00619 2.122e-5], [0.002301 7.897e-6],Ts) ;


for ct=1:N
%     ct
    u1=ones(2,1)*[ rVector(ct,:),y];
    [controls, t, x1] = lsim(CLry,u1,0:1,x);
    [controlsI, tI, x1I] = lsim(CLryI,u1,0:1,xI);

    
    x=x1(2,:);
    xI=x1I(2,:);
    
    uH=controls(2,1)+actionIFh2;
    uT=controls(2,2)+actionIFc2;
    cVector(ct,1)=max(min(uH,UB),LB);
    cVector(ct,2)=max(min(uT,UB),LB);

    
    % error
    e = rVector(ct,:) - y;
    
    % anti windup
    if UseBackCalculation
        actionIFh2 = actionIFh2 + (KiFh*e(1) + (cVector(ct,1)-uH)/tau)*Ts;
        actionIFc2 = actionIFc2 + (KiFc*e(2) + (cVector(ct,2)-uT)/tau)*Ts;
    else
%         actionIFh = actionIFh + KiFh*e(1)*Ts;
%         actionIFc = actionIFc + KiFc*e(2)*Ts;
    end
    





    
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
    uVecFh(ct+delayC) = u(1,1);
    uVecFc(ct) = u(1,2);
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

figure
plot(cVector(:,1),'--r')
hold on
plot(cVector(:,2),'--g')
title("u")
hold off


