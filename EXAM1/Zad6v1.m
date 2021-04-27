clc
clear all
close all

N=4;
Nu=2;
ny=2;
nu=2;

a= [0.4 0;
    0 -0.2];

b=[0.2 0.1;
   0.1 0.3];
  
yVector = zeros(2,2, N);
uVector1= [1;0];
uVector2= [0;1];


for i=2:N
    yVector(1, 1, i)=(a(1,1)*yVector(1, 1, i-1)+b(1,1));
    yVector(1, 2, i)=(a(1,1)*yVector(1, 2, i-1)+b(1,2));
    yVector(2, 1, i)=(a(2,2)*yVector(2, 1,i-1)+b(2,1));
    yVector(2, 2, i)=(a(2,2)*yVector(2, 2,i-1)+b(2,2));
      
end

% (N,Nu,ny,nu,steps)
[M, Mp]=DMCmatrices(yVector, N,Nu)