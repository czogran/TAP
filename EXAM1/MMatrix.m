function [M] = MMatrix(N,Nu,ny,nu,steps)
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
             Smatrix(x,y)=S((x-1)*2+y);
          end
        end

        M((i-1)*ny+1:i*ny,(j-1)*nu+1:j*nu)=Smatrix;
    end
end
end

