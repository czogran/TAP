transferModel

N=300;
Nu=50;
simLength=3000;
Ts=Tp;

nu=2;
ny=2;
psi=1;
lambda=1;

% control boundaries
duMin = -1;
duMax = 1;
uMin=0;
uMax=100;

G=transmit;

% function [M,CtAt,CtV]=MPCSmatrices(A,B,C,N,Nu)
[M,CtAt,CtV]=MPCSmatrices(A,B,C,N,Nu);
[M1,CtAt1,CtV1]=MPCSmatrices(G.A,G.B,G.C,N,Nu);

GDiscrete=c2d(G,Tp);
% GDiscrete=G;

AD=GDiscrete.A;
BD=GDiscrete.B(:,1:2);
BDAll=GDiscrete.B;
CD=GDiscrete.C;
[M,CtAt,CtV]=MPCSmatrices(AD,BD,CD,N,Nu);

psiMatrix=eye(ny*N)*psi;
lambdaMatrix=eye(nu*Nu)*lambda;

K=(M'*psiMatrix*M+lambdaMatrix)^(-1)*M'*psiMatrix;

K1=K(1:nu,:);
Ke=0;
for i=1:N
    Ke=+K1(:,(i-1)*ny+1:i*ny);
end



% initial condition
r1 = [V0, T0];
r2=[V0-2000, T0-10];
rVector=[ones(2000,2).*r1;ones(2000,2).*r2];


% vector pre allocation
uVec=zeros(simLength,6);
yVec=zeros(ny,simLength);
y=[0;0];
uFh=zeros(1,simLength);
uFc=zeros(1,simLength+delayC);
hVector=zeros(1,simLength);

inputs=zeros(6,1);

for k=3:simLength
    yZk=rVector(k,:)';
%     k
%     y
%     AD
%     yVec(:,k-1)
%     BD
%     uVec(:,k-1)
    yPrev=[(yVec(1,k-1));yVec(2,k-1)];
    yPrev2=[(yVec(1,k-2));yVec(2,k-2)];
    v = yPrev-AD*yPrev2-BDAll*uVec(k-2,:)';
  
    du=0;
   
    Vp=eye(ny);
    for i=1:N
        Kp=K1(:,(i-1)*ny+1:i*ny);
%         Ap=CtAt((i-1)*ny+1:i*ny,:);
%         Vp=CtV((i-1)*ny+1:i*ny,:);
        
%         du=du+Kp*(yZk -Ap*y-Vp*BD*uVec(:,k-1)-Vp*v);
        du=du+Kp*(yZk-C*AD^i*yPrev-C*Vp*BDAll*uVec(k-2,:)'-C*Vp*v);
        Vp=Vp+A^(i);

    end
%     du=Ke*(yZk)-K1*CtAt*y-K1*CtV*(BD*uVec(:,k-1)-v);

%         du=-K1*CtAt*y-K1*CtV*(-v);    
    
%     y_zad_act = zeros(1,N*ny);
%     for i=2:2:2*N
%        y_zad_act(1,i-1) = 1;
%        y_zad_act(1,i) = 1;
%     end
%     du=K1*(y_zad_act' - CtAt*y - CtV * [BD,GDiscrete.B(:,3:6)] * [uVec(:,k-1);Fd0;Th0;Tc0;Td0] - CtV*v);

    du(1)=max(min(du(1),duMax), duMin);
    du(2)=max(min(du(2),duMax), duMin);
    
    
    uFh(k)=max(min(uFh(k-1)+du(1),uMax),uMin);
    uFc(k+delayC)=max(min(uFc(k-1+delayC)+du(2),uMax),uMin);

    V=y(1);
    Finputs=[[uFh(k),uFc(k)],Fd];
    Tinputs=[Th,Tc,Td];
    
    uVec(k,:)=[Finputs,Tinputs];
    
    if(k*Tp<delayC)
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
    y(1)=(V);
    hVector(k)=V;
            
    T=y(2);
    kT1= dTdt(V,T,delay,Finputs,Tinputs);
    kT2= dTdt(V + Ts*V/2,T + Ts/2*kT1,delay,Finputs,Tinputs);
    kT3= dTdt(V + Ts*V/2,T + Ts/2*kT2,delay,Finputs,Tinputs);
    kT4= dTdt(V + Ts*V,T + Ts*kT3,delay,Finputs,Tinputs);       
    
    dT=Ts/6*(kT1+2*kT2+2*kT3+kT4);
    y(2)=T+ dT;
    
    yVec(1,k) = ( y(1));
    yVec(2,k+delayT) = y(2);
%     uFh(k) = u(1,1);
%     uFc(k+delayC) = u(1,2);
end

figure
plot(yVec(1,:))
hold on
plot(rVector(:,1))
title("V")
hold off

figure
plot(yVec(2,:))
hold on
plot(rVector(:,2))
title("T")
hold off

figure
plot(uFh)
hold on
plot(uFc)
legend("uFh", "uFc")
title("u")
hold off