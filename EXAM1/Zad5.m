clear all;
s(1,1,:)=[0, 0.4, 0.8, 1, 1.1, 1.15, 1.175];



dim=size(s);
ny=dim(1);
nu=dim(2);
D = dim(3);
temp = zeros(ny,nu);
S = zeros(ny,nu,D);
for i = 1:D
    for nyy = 1:ny
        for nuu = 1:nu
            temp(nyy,nuu) = s(nyy,nuu,i);
        end
    end
    S(:,:,i)=temp; 
end
%parametry
N=3;
Nu=2;
[M, MP] = DMCmatrices(S,N,Nu);

dim_M = size(M);
%parametry
psi = repmat([1],dim_M(1)/ny,1);
lambda = repmat([0.1],dim_M(2)/nu,1);

Psi = diag(psi);
Lambda = diag(lambda);
H=2*(M'*Psi*M+Lambda);

%Macierz K
K=(M'*Psi*M+Lambda)^(-1)*M';
K1=K(1,:);
% wzmocnienie uchybu
Ke=sum(K1);
y = 9;
yzad=20;
% du od uchybu

% wektor wzmocnien poprzednich sterowañ
Ku=K1*MP;