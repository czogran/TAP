clc;
close all;
clear all;

N=3;
Nu=2;
psi=1;
lambda=0.5;

ny=2;
nu=2;

psiMatrix=eye(ny*N)*psi;
lambdaMatrix=eye(Nu*nu)*lambda;

A=[0.4 0.2;
   0.1 0.5];

B =[1 0;
    0 1];

C=[1  0;
   0  1];

[M,CtAt,CtV]=MPCSmatrices(A,B,C,N,Nu);

K=(M'*psiMatrix*M+lambdaMatrix)^(-1)*M'*psiMatrix;
K1=K(1:nu,:);

% TO JEST ZLE!!!!
Ke=zeros(nu,ny);
for i=1:nu
    for j=1: ny
        sum(K1(i,(j-1)*N+1:j*N))
        Ke(i,j)=sum(K1(i,(j-1)*N+1:j*N));
    end
end

s=0;
for i=1:N
   s=s+ K(1:2,(i-1)*ny+1:i*ny)*C*A^i
end


itemsX=K1*CtAt;
itemsU=K1*CtV*B;
itemsV=K1*CtV;

seed=[0.1 0.15];
Kobs=place(A',C',seed);
