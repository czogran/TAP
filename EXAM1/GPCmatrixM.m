function M=GPCmatrixM(A,B,N,Nu)
% macierz dynamiczna M dla modelu MIMO y(k)=-A*y(k)+B*u(k-1),
% A - macierz diagonalna wym. ny x ny x nA; B - macierz wym. ny x nu x nB+1,
% N - dla horyzont predykcji; Nu - horyzont sterowania
ny=size(B,1); nu=size(B,2); nB=size(B,3)-1; nA=size(A,3);
% Współczynniki odpowiedzi skokowej macierzowej na horyzoncie N:
S=zeros(ny,nu,N);
for ks=1:N
    for m=1:ny
        for j=1:nu
            sx=0; 
            for i=1:min(ks-1,nA)
                sx=sx-A(m,m,i)*S(m,j,ks-i); 
            end
            for i=1:min(ks,nB+1)
                sx=sx+B(m,j,i); 
            end
            S(m,j,ks)=sx;
        end
    end
end

M=zeros(N*ny,Nu*nu); % Macierz dynamiczna
for i=1:N, M((i-1)*ny+1:(i-1)*ny+ny,1:nu)=S(:,:,i); end
for i=2:Nu
    M(:,(i-1)*nu+1:(i-1)*nu+nu)=[zeros((i-1)*ny,nu); M(1:(N-i+1)*ny,1:nu)];
end