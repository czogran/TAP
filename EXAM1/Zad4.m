clear all;
s(1,1,:)=[1, 2, 2.5, 2.6, 2.6];
s(1,2,:)=[0.5, 1, 1.3, 1.5, 1.5];
s(2,1,:)=[0.6, 1.2, 1.5, 1.6, 1.6];
s(2,2,:)=[0.5, 1.2, 1.6, 1.8, 1.8];

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
psi = repmat([1; 1],dim_M(1)/ny,1);
lambda = repmat([0.1; 0.1],dim_M(2)/nu,1);

Psi = diag(psi);
Lambda = diag(lambda);
H=2*(M'*Psi*M+Lambda)
