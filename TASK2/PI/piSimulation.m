% PI regulator
init


% CONFIG
% decoupler
useDecoupler =true;
% disturbance
disturbance = true;
% use linear object model
useLinearModel = false;


UseBackCalculation = true; % enable/disable anti-windup (back calculation)
KpFh = 0.9; KiFh = 0.01; % PI controller gains (parallel)
KpFc = -0.2; KiFc = -0.002; % PI controller gains (parallel)

Ts = 1; % controller sample time
tau = 1; % reset time constant
UB = 300; LB = 0; % saturation limits

% closed-loop simulation (200 steps)
N = 18000;

% initial condition
interval=N/6;
r1 = ones(interval,2).*[h0, T0];
r2=ones(interval,2).*[h0-10, T0];
r3=ones(interval,2).*[h0+10, T0];
r4=ones(interval,2).*[h0,T0-10];
r5=ones(interval,2).*[h0,T0+10];
r6=ones(interval,2).*[h0+10,T0-10];
rVector=[r1;r2;r3;r4;r5;r6];

% Decoupler
yD21=0.002301;
y1D21=7.897e-6;
uD21=0.00619;
u1D21=2.122e-5;

if disturbance
    
else
    Fdvector= ones(N,1)*Fd0;
    Tdvector= ones(N,1)*Td0;
end

if useLinearModel & useDecoupler
elseif useLinearModel &  ~useDecoupler
elseif ~useLinearModel &  ~useDecoupler
    piNoLinearNoDecouplerBody
elseif ~useLinearModel &  useDecoupler
    piNoLinearDecouplerBody
       
end
