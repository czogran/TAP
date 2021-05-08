transferModelNoDelays;

UseBackCalculation = true; % enable/disable anti-windup (back calculation)
G=tf(transmit);
KpFh = 10; KiFh = 0.1; % PI controller gains (parallel)
KpFc = -0.5; KiFc = -0.01; % PI controller gains (parallel)


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
rVector=[ones(N-2000,2).*r1;ones(2000,2).*r2];

% vector pre allocation
y = [0,0];

offset=max(delayC,delayT);
yVec = zeros(N+offset,2);
uVecFh= zeros(N+offset,1);
uVecFc= zeros(N+offset,2);

% Decoupler
yD21=0.002301;
y1D21=7.897e-6;
uD21=0.00619;
u1D21=2.122e-5;
D21=tf([0.00619 2.122e-5], [0.002301 7.897e-6],Ts) ;

uVecD21= zeros(N+offset,2);
yVecD21= zeros(N+offset,2);
dOutputs=zeros(N+offset,1);

for ct=1:N
    % error
    e = rVector(ct,:) - y;
    
    % control action
    actionPFh = KpFh*e(1);
    actionPFc = KpFc*e(2);

    uFh = actionPFh + actionIFh;
    uFc = actionPFc + actionIFc;

    % saturation control action
    uFhSat = max(min(uFh,UB),LB);
    uFcSat = max(min(uFc,UB),LB);
    
    yD12=uFcSat*(-1);
    uFhSatDHelp=uFhSat+yD12;
    uFhSatD=max(min(uFhSatDHelp,UB),LB);
    
    yD21=(uD21*uFhSat+u1D21*uVecD21(max(ct-1,1))-y1D21*yVecD21(max(ct-1,1)))/yD21;
%     yD21=uFhSat*2.69;

    uFcSatDHelp=uFcSat+yD21;
    uFcSatD=max(min(uFcSatDHelp,UB),LB);
    
    if(ct<130)
        ct
        uFhSat
        yD12
        uFhSatDHelp
        uFhSatD
%         (uD21*uFhSat)
%         u1D21*uVecD21(max(ct-1,1))
%         y1D21*yVecD21(max(ct-1,1))
        
    end
    
    uVecD21(ct)=uFhSat;
    yVecD21(ct)=yD21;
    
    u=[uFhSatD, uFcSatD];
    
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