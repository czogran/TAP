function [Mp] = MpMatrix(N,D,ny,nu,steps)
Mp= zeros(N*ny,(D-1)*nu);
index=0;
for i =1:N
    for j =1:D-1
       for x=1:ny
           for y=1:nu 
                if i+j<=D
                   Mp(ny*(i-1)+x,nu*(j-1)+y)=steps(2*(x-1)+y,j+i)-steps(2*(x-1)+y,j);
                else
                   Mp(ny*(i-1)+x,nu*(j-1)+y)=steps(2*(x-1)+y,D)-steps(2*(x-1)+y,j);
                end
               index=index+1;
%                Mp(ny*(i-1)+x,nu*(j-1)+y)=steps(2*(x-1)+y,i+1)-steps(2*(x-1)+y,j);      
           end
       end
    end
end
end

