clear all;

N=5;
Nu=2;
ny=1;
nu=1;
psi=1;
lambda=0.25;

yZad=20;
y=9;

steps = [ 0 2 3 3.5 4];
stepsZ= [1 1.5 1.75 2 2];

D=length(steps);

u=[4 3 1 0 -4];
z=[5 1 1 3 1];

du= circshift(u,1)- u;
du=du(2:D)

dz=circshift(z,1)- z;
dz=[dz(2:D),0]

% MMatrix(N,Nu,ny,nu,steps)
M=MMatrix(N,Nu,ny,nu,steps);
% MpMatrix(N,D,ny,nu,steps)
Mp=MpMatrix(N,D,ny,nu,steps)
MZp=MpMatrix(N,D,ny,nu,stepsZ)
MZp=[stepsZ',MZp]


psiMatrix= eye(N)*psi;
lambdaMatrix =  eye(Nu)*lambda;

K= (M'*psiMatrix*M+lambdaMatrix)^(-1)*M'*psiMatrix
Ke=sum(K(1,:))

K1=K(1,:)
K1*Mp


duu = Ke*(yZad-y)-K1*Mp*du'

duz = duu  - K1*MZp*dz'