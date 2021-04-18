
N=3;
Nu=2;

ny=2;
nu=2;

psiValues=[1,1];
lambdaValues=[0.1,0.1];


steps = [ 1 2 2.5 2.6 2.6;
    0.5 1 1.3 1.5 1.5;
    0.6 1.2 1.5 1.6 1.6;
    0.5 1.2 1.6 1.8 1.8];

psiVector= zeros(N*ny,1);
for i=1:ny
   psiVector(((i-1)*N+1):i*N)=psiValues(i);
end

lambdaVector= zeros(Nu*nu,1);
for i=1:nu
   lambdaVector(((i-1)*Nu+1):i*Nu)=lambdaValues(i);
end

psi=diag(psiVector)

lambda=diag(lambdaVector)

M= zeros(N*ny,Nu*nu);

for i =1:N
    for j =1:Nu
        if(i<j)
            continue;
        end
        index=i-j+1;
        S=steps(:,index);
        Smatrix= zeros(nu,ny);
        for x=1:ny
          for y=1:nu
             Smatrix(x,y)=S((x-1)*2+y)
          end
        end
        
        M((i-1)*ny+1:i*ny,(j-1)*nu+1:j*nu)=Smatrix;
    end
end



H = 2*(M'*psi*M+lambda)
