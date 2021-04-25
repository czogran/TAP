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
  
yVector = zeros(4,N);
uVector1= [1;0];
uVector2= [0;1];


for i=2:N
    s11=(a*yVector(1,i-1)+b*uVector1);
    s12=(a*yVector(1,i-1)+b*uVector2);
    s21=(a*yVector(1,i-1)+b*uVector1);
    s22=(a*yVector(1,i-1)+b*uVector2);
       
    yVector(1,i) = s11(1);
    yVector(2,i) = s12(1);
    yVector(3,i) = s21(1);
    yVector(4,i) = s22(2);

end

% (N,Nu,ny,nu,steps)
M=MMatrix(N,Nu,ny,nu,yVector)