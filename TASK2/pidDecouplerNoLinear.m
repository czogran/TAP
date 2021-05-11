% PI regulator
init

UseBackCalculation = true; % enable/disable anti-windup (back calculation)
KpFh = 0.04; KiFh = 0.0004; % PI controller gains (parallel)
KpFc = 0; KiFc = -0.01; % PI controller gains (parallel)

Ts = Tp; % controller sample time
tau = 1; % reset time constant
UB = 1500; LB = 0; % saturation limits

% closed-loop simulation (200 steps)
N = 18000;

% initial condition
interval=N/6;
r1 = ones(interval,2).*[h0, T0];
r2=ones(interval,2).*[h0-10, T0];
r3=ones(interval,2).*[h0+10, T0];
r4=ones(interval,2).*[h0,T0-10];
r5=ones(interval,2).*[h0,T0+5];
r6=ones(interval,2).*[h0+10,T0-5];
rVector=[r1;r1;r6;r5;r5;r5];

% vector and values pre allocation
y = [0,0];

offset=max(delayC,delayT);
yVec = zeros(N+offset,2);
uVecFh= zeros(N+offset,1);
uVecFc= zeros(N+offset,2);
uVecFcin= zeros(N+offset,2);

uVecD21= zeros(N+offset,1);
yVecD21= zeros(N+offset,2);

uVecD12= zeros(N+offset,1);
yVecD12= zeros(N+offset,2);

actionIFh = 0;
actionIFc = 0;

errorH=0;
errorT=0;

% Decoupler
yD21Const=0.002301;
y1D21Const=7.897e-6;
uD21Const=0.00619;
u1D21Const=2.122e-5;
for ct=1:N
    % error
    e = rVector(ct,:) - y;
    errorH=errorH+e(1)^2;
    errorT=errorT+e(2)^2;
    
    % control action
    actionPFh = KpFh*e(1);
    actionPFc = KpFc*e(2);

    uFh = actionPFh + actionIFh;
    uFc = actionPFc + actionIFc;

    % saturation control action
    uFhSat = max(min(uFh,UB),-UB);
    uFcSat = max(min(uFc,UB),-UB);
    
    
    % anti windup
    if UseBackCalculation
        actionIFh = actionIFh + (KiFh*e(1) + (uFhSat-uFh)/tau)*Ts;
        actionIFc = actionIFc + (KiFc*e(2) + (uFcSat-uFc)/tau)*Ts;
    else
        actionIFh = actionIFh + KiFh*e(1)*Ts;
        actionIFc = actionIFc + KiFc*e(2)*Ts;
    end
    
    % decoupler     
    yD12=uFcSat*(-1);
    uFhSatDHelp=uFhSat+yD12;
    uFhSatD=max(min(uFhSatDHelp,UB),LB);
    
    yD21=(uD21Const*uFhSat+u1D21Const*uVecD21(max(ct-1,1))-y1D21Const*yVecD21(max(ct-1,1)))/yD21Const;
%     yD21=uFhSat*2.69;
    uFcSatDHelp=uFcSat+yD21;
    uFcSatD=max(min(uFcSatDHelp,UB),LB);

    uVecFh(ct) = uFhSatD;
    uVecFc(ct+delayC/Ts) = uFcSatD;
    uVecFcin(ct) = uFcSatD;
    
    uVecD21(ct)= uFhSat;
    yVecD21(ct)=yD21;
    
    uVecD12(ct)=uFcSat;
    yVecD12(ct)=yD12;
 
    u=[uVecFh(ct), uVecFc(ct)];
    
    V=volume(y(1));
    Finputs=[u,Fd];
    Tinputs=[Th,Tc,Td];
    
    delay=1;
 
    kV1= dVdt(heightFromVolume(V), delay, Finputs);
    kV2= dVdt(heightFromVolume(V + Ts/2*kV1),delay,Finputs);
    kV3= dVdt(heightFromVolume(V + Ts/2*kV2),delay,Finputs);
    kV4= dVdt(heightFromVolume(V + Ts*kV3),delay,Finputs);
    
    dV=Ts/6*(kV1+2*kV2+2*kV3+kV4);
    V=V+dV;
    % output h     
    y(1)=heightFromVolume(V);
            
    T=y(2);
    kT1= dTdt(V,T,delay,Finputs,Tinputs);
    kT2= dTdt(V + Ts*V/2,T + Ts/2*kT1,delay,Finputs,Tinputs);
    kT3= dTdt(V + Ts*V/2,T + Ts/2*kT2,delay,Finputs,Tinputs);
    kT4= dTdt(V + Ts*V,T + Ts*kT3,delay,Finputs,Tinputs);       
    
    dT=Ts/6*(kT1+2*kT2+2*kT3+kT4);
    %  output T   
    y(2)=T+ dT;
    
    % outputVector     
    yVec(ct,1) = y(1);
    yVec(ct+delayT,2) = y(2);
end

% remove offset
yVec = yVec(1:N,:);
uVecFh= uVecFh(1:N);
uVecFcin= uVecFcin(1:N);
uVecFc= uVecFc(1:N);


clf;
figure(1);
plot(1:N,yVec(:,1),'r');
hold on
plot(rVector(:,1),'--b');
xlabel('Time'); ylabel('Signal');
title("Regulator PI z odsprzęganiem"+newline+"h[cm]" );
xlabel("t[s]")
ylabel("h[cm]")
legend("h[cm]", "trajektoria zadana", 'Location','best')
hold off

figure(2);
plot(1:N,yVec(:,2),'.r');
hold on
plot(rVector(:,2),'--b');
title("Regulator PI z odsprzęganiem"+newline+"T[\circC]");
xlabel('t[s]'); ylabel('T[\circC]');
legend('temperatura [\circC]','trajektoria zadana', 'Location','best')
hold off

figure
plot(uVecFh,'--r')
hold on
plot(uVecFcin,'.g')
title("u")
title("Regulator PI z odsprzęganiem"+newline+"sterowanie wejściowe obiektu"+newline+"sterowanie u");
legend("Fh","Fcin", 'Location','best')
xlabel("t[s]")
ylabel("u[$\frac{cm^3}{s}$]",'Interpreter','latex')
hold off


figure
plot(uVecD21,'--r')
hold on
plot(uVecD12,'--g')
title("u")
title("Regulator PI z odsprzęganiem"+newline+"sterowanie wyjściowe regulatorów"+newline+"sterowanie u");
legend("PI 1","PI 2", 'Location','best')
xlabel("t[s]")
ylabel("u[$\frac{cm^3}{s}$]",'Interpreter','latex')
hold off
