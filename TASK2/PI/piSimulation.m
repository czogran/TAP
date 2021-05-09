% PI regulator
init

% CONFIG
% decoupler
useDecoupler =true;
% disturbance
disturbance = true;
% use linear object model
useLinearModel = false;
% dteremines for 0 simulations start from zero, if 1 it starts from work point
initFactor=1;

UseBackCalculation = true; % enable/disable anti-windup (back calculation)

Ts = 1; % controller sample time
tau = 1; % reset time constant
UB = 300; LB = 0; % saturation limits

% closed-loop simulation (200 steps)
N = 18000;

% initial condition
interval=N/3;
r1 = ones(interval,2).*[h0, T0];
r2=ones(interval,2).*[h0-10, T0];
r3=ones(interval,2).*[h0+10, T0];
r4=ones(interval,2).*[h0,T0-10];
r5=ones(interval,2).*[h0,T0+10];
r6=ones(interval,2).*[h0+10,T0-5];
rVector=[r1;r5;r6];


% path for saving
fileName="step-responses";
overLeafFilePath="img/PI/";
path="../img/PI/";

if useDecoupler
    overLeafFilePath=overLeafFilePath+"decoupler/";
    path=path+"decoupler/";
    
    % PI params with decoupler
    KpFh = 0.04; KiFh = 0.0004; % PI controller gains (parallel)
    KpFc = 0; KiFc = -0.01; % PI controller gains (parallel)

    % Decoupler params
    yD21Const=0.002301;
    y1D21Const=7.897e-6;
    uD21Const=0.00619;
    u1D21Const=2.122e-5;
else
    overLeafFilePath=overLeafFilePath+"noDecoupler/";
    path=path+"noDecoupler/";
    
    % PI without decoupler      
    KpFh = 0.9; KiFh = 0.01; % PI controller gains (parallel)
    KpFc = -0.2; KiFc = -0.002; % PI controller gains (parallel)
end

if disturbance
    
else
    Fdvector= ones(N,1)*Fd0;
    Tdvector= ones(N,1)*Td0;
end

if useLinearModel 
    T0inputs=[Th0,Tc0,Td0];
    F0inputs=[Fh0,Fcin0,Fd0];
end
    
if useDecoupler
    index=1;
    piDecouplerBody;
else
    index=1;
    piNoDecouplerBody ;      
end
