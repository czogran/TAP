clear all;
s(1,1,:)=[0, 2, 3, 3.5, 4];
z(1,1,:)=[1, 1.5, 1.75, 2, 2];
zz=[1, 1.5, 1.75, 2, 2];

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
            temp_z(nyy,nuu) = z(nyy,nuu,i);
        end
    end
    S(:,:,i)=temp; 
    Z(:,:,i)=temp_z;
end
%parametry
N=5;
Nu=2;
[M, MP] = DMCmatrices(S,N,Nu);
[MZ, MZP] =DMCmatrices(Z,N,Nu);
MZP = [zz',MZP];
dim_M = size(M);
%parametry
psi = repmat([1],dim_M(1)/ny,1);
lambda = repmat([0.25],dim_M(2)/nu,1);

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
duy = Ke*(yzad-y);
% wektor wzmocnien poprzednich sterowañ
Ku=K1*MP;
du = [1, 2, 1, 4];
% du od poprzednich sterowañ
duu=sum(Ku.*du);
% wektor wzmocnieñ zak³óceñ
Kuz=K1*MZP;
dz=[4, 0, -2, 2, 0];
% du od poprzednich zak³óceñ
duz=sum(Kuz.*dz);

%prawo regulacji:
du = duy - duu - duz