transferModel

N=400;
Nu=50;
simLength=20000;
Tp;

nu=2;
ny=2;
psi=1;
lambda=1;

% control boundaries
duMin = -0.05;
duMax = 0.05;
uMin=0;
uMax=300;

G=transmit;
GDiscrete=c2d(G,Tp);
% GDiscrete=G;

AD=GDiscrete.A;
BD=GDiscrete.B(:,1:2);
BDAll=GDiscrete.B;
CD=GDiscrete.C;
[M,CtAt,CtV]=MPCSmatrices(AD,BD,CD,N,Nu);

psiMatrix=eye(ny*N)*psi;
lambdaMatrix=eye(Nu*nu)*lambda;

K=(M'*psiMatrix*M+lambdaMatrix)^(-1)*M'*psiMatrix;

K1=K(1:nu,:);

% initial condition
r1 = [1000, -5];
r2=[V0-2000, T0-10];
rVector=[zeros(simLength/2,2);ones(simLength/2,2).*r1];

x0=[V0;T0];

% vector pre allocation
uVecIn=zeros(simLength,6);
uVec=zeros(simLength,6);

yVec=zeros(ny,simLength);
x0=[V0;T0];
xVec=ones(ny,simLength).*x0;

y=[V0;T0];
uFh=ones(1,simLength)*Fh0;
uFcin=ones(1,simLength+delayC)*Fc0;
uFc=ones(1,simLength+delayC)*Fc0;

hVector=zeros(1,simLength);

inputs=zeros(6,1);

du=zeros(1,6);

for k=3:simLength
    yZk=rVector(k,:);
    
    uPrev=uVecIn(k-1,:)'-[Fh0;Fc0;Fd0;Th0;Tc0;Td0];
%     uPrev=zeros(6,1);

    x=xVec(:,k)-x0;
    xPrev=xVec(:,k-1)-x0;
    v = x-AD*xPrev-BDAll*uPrev;
%     v=[0;0];
    
%     dU=K1*(yZk'-CtAt*yPrev-CtV*BDAll*uVec(k-2,:)'-CtV*v);
    
    dU=K1*( repmat(yZk',N,1)-CtAt*x-CtV*BDAll*uPrev-CtV*v);

    du(1)=max(min(dU(1),duMax), duMin);
    du(2)=max(min(dU(2),duMax), duMin);
    
    uFh(k)=max(min(uFh(k-1)+du(1),uMax),uMin);
    uFcin(k)=max(min(uFcin(k-1)+du(2),uMax),uMin);
    uFc(k+delayC)=max(min(uFc(k-1+delayC)+du(2),uMax),uMin);

    V=y(1);
    
%     TODO jak jest Fcin to smiga
    Finputs=[[uFh(k),uFc(k)],Fd];
    Tinputs=[Th,Tc,Td];
    
    uVecIn(k,:)=[Finputs,Tinputs];
    uVec(k,:)=[[[uFh(k),uFc(k)],Fd],Tinputs];

    delay=1;
    
    kV1= dVdt(heightFromVolume(V), delay, Finputs);
    kV2= dVdt(heightFromVolume(V + Tp/2*kV1),delay,Finputs);
    kV3= dVdt(heightFromVolume(V + Tp/2*kV2),delay,Finputs);
    kV4= dVdt(heightFromVolume(V + Tp*kV3),delay,Finputs);
    
    dV=Tp/6*(kV1+2*kV2+2*kV3+kV4);
    V=V+dV;
    y(1)=(V);
    hVector(k)=V;
            
    T=y(2);
    kT1= dTdt(V,T,delay,Finputs,Tinputs);
    kT2= dTdt(V + Tp*V/2,T + Tp/2*kT1,delay,Finputs,Tinputs);
    kT3= dTdt(V + Tp*V/2,T + Tp/2*kT2,delay,Finputs,Tinputs);
    kT4= dTdt(V + Tp*V,T + Tp*kT3,delay,Finputs,Tinputs);       
    
    dT=Tp/6*(kT1+2*kT2+2*kT3+kT4);
    y(2)=T+ dT;
    
    yVec(1,k) = ( y(1));
    yVec(2,k+delayT) = y(2);
    
    xVec(:,k+1)=[V;T];
end

figure
plot(yVec(1,:))
hold on
plot(rVector(:,1)+V0)
title("V")
hold off

figure
plot(yVec(2,:))
hold on
plot(rVector(:,2)+T0)
title("T")
hold off

figure
plot(uFh)
hold on
plot(uFc)
legend("uFh", "uFc")
title("u")
hold off