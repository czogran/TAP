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
r6=ones(interval,2).*[h0+10,T0-10];
r7=ones(interval,2).*[h0-10,T0+5];

rArray={[r1;r6;r7],[r1;r2;r3], [r1;r4;r5]};
rWorkPoint=[r1;r1;r1];
% rVector=[r1;r5;r6];


% path for saving
fileName="pi";
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
    overLeafFilePath=overLeafFilePath+"disturbance/";
    path=path+"disturbance/";
    
    rFd1 = ones(interval,1)*Fd0;
    rFd2=ones(interval,1)*(Fd0+10);
    rFd3=ones(interval,1)*(Fd0-10);

    rTd1 = ones(interval,1)*(Td0);
    rTd2=ones(interval,1)*(Td0+10);
    rTd3=ones(interval,1)*(Td0-10);
    
    % arrays with disturbance, must be same length     
    rFdArray={[rFd1 ;rFd2 ;rFd3] [rFd1; rFd1; rFd1]};
    rTdArray={[rTd1 ;rTd1 ;rTd1] [rTd1; rTd2 ;rTd3]};
else
    overLeafFilePath=overLeafFilePath+"noDisturbance/";
    path=path+"noDisturbance/";
    
    FdVector= ones(N,1)*Fd0;
    TdVector= ones(N,1)*Td0;
end

% setup work point vectors for linear model
if useLinearModel 
     fileName=fileName+"LinearModel";
    
    T0inputs=[Th0,Tc0,Td0];
    F0inputs=[Fh0,Fc0,Fd0];
else
     fileName=fileName+"NoLinearModel";
end
   

if(initFactor==0)
    fileName=fileName+"ZeroStart"
    
    index = 0
    rVector=rWorkPoint;
    if useDecoupler
        piDecouplerBody;
    else
        piNoDecouplerBody
    end
elseif disturbance
     iterator=1;
    for i =1: length(rFdArray)
        iterator=iterator+1;
        index= iterator+("Dist"+disturbance)+("Lin"+useLinearModel);
        
        FdVector = rFdArray{i};
        TdVector = rTdArray{i};
        rVector=rWorkPoint;
        if useDecoupler
            piDecouplerBody;
        else
            piNoDecouplerBody;
        end
    end
else
    iterator=1;
    for i =1: length(rArray)
        index= iterator+("Lin"+useLinearModel);
        iterator=iterator+1;

        rVector = rArray{i};
        if useDecoupler
            piDecouplerBody;
        else
            piNoDecouplerBody;
        end
    end
end

