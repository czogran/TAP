clear all;
close all;

transferModel;
sys_con=transmit;
Ts=Tp;
Tsim = 3000;
lambda_coefficient = 1;
C_coefficient=0.7;


sys_dyskr = c2d(sys_con, Ts);
% sys_dyskr=transmit;

N = 300;
Nu = 50;

% restrictions on control signal.
uMin=[-5,-5];
deltaUMin=[-0.5,-0.5];
deltaUMax=[0.5,0.5];
uMax=[5, 5];

% restrictions on outputs.
xMin=[-inf,-inf];
xMax=[inf,inf];

teta = obsv(sys_dyskr);
A_dyskr =sys_dyskr.A;
B1_dyskr= sys_dyskr.B(:,1:2);
B2_dyskr= sys_dyskr.B(:,3:4);
C_dyskr=sys_dyskr.C;
D_dyskr=sys_dyskr.D;
[M,CtAt,CtV]=MPCSmatrices(A_dyskr, B1_dyskr, C_dyskr, N, Nu);

[nx, nu]=size(B1_dyskr); ny=size(C_dyskr,1);

psi=eye(N*ny);
lambda=eye(Nu*nu)*lambda_coefficient;
K=(M'*psi*M+lambda)^-1*M'*psi;
K1=K(1:nu,:);

% if rank(teta) == nx
%     
% else
%     Kobsp = place(A_dyskr', C_dyskr', zeros(1, nx));
%     Kobsb = A^-1 * Kobsp;
% end

%%Simulation
sim_len = length(1:Ts:Tsim);
u=zeros(nu,sim_len);
y=zeros(ny,sim_len);
x=zeros(nx,sim_len);
v=zeros(nx,sim_len);
y_zad=[zeros(ny,1000),ones(ny,sim_len-1000)];
%y_zad=[zeros(ny,1000),zeros(ny,sim_len-1000)];
FD = ones(1,sim_len)*Fd0;
TD = ones(1,sim_len)*Td0;

h = ones(1,sim_len)*h0;
T = ones(1,sim_len)*Tp;
Tout = ones(1,sim_len)*Tp;


for k = 200 : sim_len-N
    %Prepare y_zad.
    y_zad_act = zeros(1,N*ny);
    for i=2:2:2*N
       y_zad_act(1,i-1) = y_zad(1,k-1+i/2);
       y_zad_act(1,i) = y_zad(2,k-1+i/2);
    end
    
    %Error.
    v(:,k)=x(:,k)-(A_dyskr*x(:,k-1)+B1_dyskr*[u(1,k-1);u(2,k-1-floor(k-delayC/Ts))] + B2_dyskr*[FD(1,k-1)-Fd0;TD(1,k-1)-Td0]);
    delta_u = K1*(y_zad_act' - CtAt*[x(1,k-floor(delayT/Ts)); x(2,k)] - CtV * [B1_dyskr, B2_dyskr] * [u(1,k-1);u(2,k-1-floor(k-delayC/Ts));FD(1,k)-Fd0;TD(1,k)-Td0] - CtV*v(:,k));
    
    %Restrictions.
    for i=1:nu
        delta_u(i,1) = min(delta_u(i,1), deltaUMax(i));
        delta_u(i,1) = max(delta_u(i,1), deltaUMin(i));
    end
    
    u(:,k) = u(:,k-1) + delta_u;
    
    for i=1:nu
        u(i,k) = min(u(i,k), uMax(i));
        u(i,k) = max(u(i,k), uMin(i));
    end
    
    %Count object
    
  h(k+1)=h(k)+1/(2*C_coefficient*h(k))*(u(1,k)+ Fh0 +u(2,k-floor(k-delayC/Ts)) +Fc0 +FD(1,k)-a*sqrt(h(k)));
  T(k+1)=T(k)+1/(C_coefficient*(h(k))^2)*((Th0-T(k))*(u(1,k) + Fh0)+(Tc0-T(k))*(u(2,k-floor(k-delayC/Ts)) + Fc0)+(TD(1,k)-T(i))*FD(1,k));

    %Restrictions
%     h(k+1)=min(h(k+1),xMax(1));
%     h(k+1)=max(h(k+1),xMin(1));        
%  
%     T(k+1)=min(T(k+1),xMax(2));
%     T(k+1)=max(T(k+1),xMin(2)); 
    
  Tout(k)=T(k-floor(delayT/Ts));
  x(:,k+1) = [T(k+1)-Tp;h(k+1)-h0];
  y(:,k) = [Tout(k)-Tp;h(k)-h0];
    
end

figure(1);
stairs(1:sim_len-N, Tout(1:sim_len-N));

figure(2);
stairs(1:sim_len-N, h(1:sim_len-N));
