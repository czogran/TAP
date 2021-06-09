clc;
close all;
clear all;

N=2;
Nu=1;
psi=1;
lambda=0.36;

ny=1;
nu=1;

psiMatrix=eye(ny*N)*psi;
lambdaMatrix=eye(Nu*nu)*lambda;

A=[0.4 0.2;
   0 0.6];

B =[0;1];

F=[0.2; 0];

C=[4 0];

[M,CtAt,CtV]=MPCSmatrices(A,B,C,N,Nu);

K=(M'*psiMatrix*M+lambdaMatrix)^(-1)*M'*psiMatrix;
K1=K(1:nu,:);
Ke=sum(K1);
itemsX=K1*CtAt;
itemsU=K1*CtV*B;
itemsV=K1*CtV;

% b
% H = 2(MT ΨM+Λ)
H=2*(M'*psiMatrix*M+lambdaMatrix);