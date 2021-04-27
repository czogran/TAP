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
    s12=(a*yVector(2,i-1)+b*uVector2);
    b*uVector1
    a*yVector(3,i-1)
    s21=(a*yVector(3,i-1)+b*uVector1)
    s22=(a*yVector(4,i-1)+b*uVector2);
       
    yVector(1,i) = s11(1);
    yVector(2,i) = s12(1);
    yVector(3,i) = s21(4);
    yVector(4,i) = s22(4);

end

% (N,Nu,ny,nu,steps)
M=MMatrix(N,Nu,ny,nu,yVector)



yVectorB = zeros(4,N);
yVectorB(:,1)=[1 1 2 2]';
uVector= [1;2];


for i=2:N
    s11=(a*yVectorB(1,i-1)+b*uVector);
    s12=(a*yVectorB(2,i-1)+b*uVector);
    s21=(a*yVectorB(3,i-1)+b*uVector);
    s22=(a*yVectorB(4,i-1)+b*uVector);
       
    yVectorB(1,i) = s11(1);
    yVectorB(2,i) = s12(1);
    yVectorB(3,i) = s21(4);
    yVectorB(4,i) = s22(4);

end

% function M=GPCmatrixM(A,B,N,Nu)
A=[0.4 0;
    0 -0.2];
B=[0.2 0.1;
    0.1 0.3];

M2=GPCmatrixM(-A,B,N,Nu)

% wektor składowej swobodnej predykcji wyjść Y0 o długości ny*N
% dla modelu obiektu MIMO y(k)=-A*y(k)+B*u(k-1), gdzie
% A - macierz diagonalna wym. ny x ny x nA; B - macierz wym. ny x nu x nB+1,
% k - bieżąca chwila czasu
% y - wektor bieżącego i poprzednich wyjść y(1:ny,1:k)
% u - wektor poprzednich sterowań u(1:nu,1:k-1)

y=[0.8 1;0.25 2];
k=2;
u=[1;2];

ny=size(B,1); nu=size(B,2); nB=size(B,3)-1; nA=size(A,3);

dk=zeros(ny,1);
for i=1:nA, dk=dk-A(:,:,i)*y(:,max(1,k-i)); end
for i=1:nB+1, dk=dk+B(:,:,i)*u(:,max(1,k-i)); end
dy=y(:,k)-dk; %estymata zakłócenia
Y0matrix=zeros(ny,N); %(macierz pomocnicza do obliczeń predykcji )
for p=1:N
yx=zeros(ny,1);
for i=1:min(nA,p-1), yx=yx-A(:,:,i)*Y0matrix(:,p-i); end
for i=min(nA,p-1)+1:nA, yx=yx-A(:,:,i)*y(:,k+p-i); end
for i=0:min(nB,p), yx=yx+B(:,:,i+1)*u(:,k-1); end
for i=min(nB,p)+1:nB, yx=yx+B(:,:,i+1)*u(:,max(1,k-1+p-i)); end
Y0matrix(:,p)=yx+dy;
end
for i=1:N, Y0((i-1)*ny+1:(i-1)*ny+ny)=Y0matrix(:,i); end

w=ss(A,B,ones(2,2),zeros(2,2));
t=20;
u1=[ones(t,1)*1,ones(t,1)*2];
lsim(w,u1,1:t)