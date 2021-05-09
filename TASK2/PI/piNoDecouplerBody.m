% PI regulator body with no linear object  and without decoupler

% vector and values pre allocation
y = [0,0];

offset=max(delayC,delayT);
yVec = ones(N+offset,2).*[V0, T0].*initFactor;
uVecFh= ones(N+offset,1)*Fh0*initFactor;
uVecFc= ones(N+offset,1)*Fc0*initFactor;
uVecFcin= ones(N+offset,1)*Fc0*initFactor;

actionIFh = 0;
actionIFc = 0;

errroH=0;
errocT=0;


for ct=1:N
    % setup disturbance     
    FdInput = FdVector(ct);
    TdInput = TdVector(ct);

    % error
    e = rVector(ct,:) - y;
    errorH=e(1)^2;
    errorT=e(2)^2;
    
    % control action
    actionPFh = KpFh*e(1);
    actionPFc = KpFc*e(2);

    uFh = actionPFh + actionIFh;
    uFc = actionPFc + actionIFc;

    % saturation control action
    uFhSat = max(min(uFh,UB),LB);
    uFcSat = max(min(uFc,UB),LB);

    uVecFh(ct) = uFhSat;
    uVecFc(ct+delayC/Ts) = uFcSat;
    uVecFcin(ct) = uFcSat;
    u=[uVecFh(ct), uVecFc(ct)];
    
    % anti windup
    if UseBackCalculation
        actionIFh = actionIFh + (KiFh*e(1) + (uFhSat-uFh)/tau)*Ts;
        actionIFc = actionIFc + (KiFc*e(2) + (uFcSat-uFc)/tau)*Ts;
    else
        actionIFh = actionIFh + KiFh*e(1)*Ts;
        actionIFc = actionIFc + KiFc*e(2)*Ts;
    end
    
    V=volume(y(1));
    Finputs=[u,FdInput];
    Tinputs=[Th,Tc,TdInput];
    
    delay=1;
 
    if useLinearModel
        linearVRk4;
    else
        noLinearVRk4;
    end
    
    dV=Ts/6*(kV1+2*kV2+2*kV3+kV4);
    V=V+dV;
    % output h     
    y(1)=heightFromVolume(V);
            
    T=y(2);
    
    if(useLinearModel)
        linearTRk4;
    else
        noLinearTRk4;
    end
    
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
uVecFcin= uVecFc(1:N);

heightFigure=figure;
plot(1:N,yVec(:,1),'r');
hold on
plot(rVector(:,1),'--b');
xlabel('Time'); ylabel('Signal');
title("Regulator PI bez odsprzęgania"+newline+"h[cm]" );
xlabel("t[s]")
ylabel("h[cm]")
legend("h[cm]", "trajektoria zadana", 'Location','best')
hold off

tempFigure=figure;
plot(1:N,yVec(:,2),'.r');
hold on
plot(rVector(:,2),'--b');
title("Regulator PI bez odsprzęgania"+newline+"T[\circC]");
xlabel('t[s]'); ylabel('T[\circC]');
legend('temperatura [\circC]','trajektoria zadana', 'Location','best')
hold off

controlPIFigure=figure
plot(uVecFh,'r')
hold on
plot(uVecFcin,'g')
title("u")
title("Regulator PI bez odsprzęgania"+newline+"sterowanie u");
legend("Fh","Fcin", 'Location','best')
xlabel("t[s]")
ylabel("u[$\frac{cm^3}{s}$]",'Interpreter','latex')
hold off

caption="Wykresy dla regulatora PI bez odsprzegania.";
label="fig:PINodDecoupler"+index;

heightName="PINoDecouplerH"+index;
tempName="PINoDecouplerT"+index;
controlPIName="PINoDecouplerControl"+index;


if ~disturbance
saveFiguresInColumn([heightFigure,tempFigure,controlPIFigure], path,[heightName,tempName,controlPIName],fileName,overLeafFilePath,caption,label);
else
    disturbanceFigure=figure
    subplot(2,1,1)
    plot(FdVector)
    title("Zakłócenie Fd")
    xlabel("t[s]")
    ylabel("Fd[$\frac{cm^3}{s}$]",'Interpreter','latex')
    subplot(2,1,2)
    plot(TdVector)
    title("Zakłócenie Td")
    xlabel("t[s]")
    ylabel("Td[\circC]")
    hold off
    sgtitle("Wartość zakłoceń")
   
    hold off


    caption="Wykresy dla regulatora PI bez odsprzegania dla różnych wartości zakłóceń";
    label="fig:PIDecoupler"+index;  
    
    controlDistubanceName="PIDecouplerDisturbance"+index;
    
    saveFiguresInColumn([heightFigure,tempFigure,controlPIFigure,disturbanceFigure], path,[heightName,tempName,controlPIName,controlDistubanceName],fileName,overLeafFilePath,caption,label);
end
