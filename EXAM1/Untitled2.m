steps = [ 1 2 2.5 2.6 2.6;
    0.5 1 1.3 1.5 1.5;
    0.6 1.2 1.5 1.6 1.6;
    0.5 1.2 1.6 1.8 1.8];


N=3;
Nu=2;

ny=2;
nu=2;
% C(D,N,Nu,lambda,psi,ny,nu,response)
% dmc= DMC(6,N,Nu,1,1,ny,nu,steps)

M=MMatrix(N,Nu,ny,nu,steps);
M
Mp1=MpMatrix(N,5,ny,nu,steps);
Mp1
% [ny,nu,D]=size(S); % D - liczba dyskretnych chwil czasu odpowiedzi skokowej
steps = [];
steps(1,1,:)=[1 2 2.5 2.6 2.6];
steps(1,2,:)=[0.5 1 1.3 1.5 1.5];

steps(2,1,:)=[0.6 1.2 1.5 1.6 1.6];
steps(2,2,:)=[0.5  1.2 1.6 1.8 1.8];

% DMCmatrices(S,N,Nu)
[M Mp] =DMCmatrices(steps, N, Nu);
Mp

M
