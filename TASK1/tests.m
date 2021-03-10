
clc
clear
close all

Tp=0.5;
Y=[0];
y=Y(1);
xVector=0:Tp:5;
for i=2:length(xVector)
    x=xVector(i-1);
    k1= fun(x,y);
    k2= fun(x+Tp/2,y+Tp/2*k1);
    k3= fun(x +Tp/2,y + Tp/2*k2);
    k4= fun(x+Tp,y+Tp*k3);

    dy=Tp/6*(k1+2*k2+2*k3+k4);
    Y(i)=Y(i-1)+dy;
    y=Y(i);
end
plot(xVector,Y)


function f = fun(x,y)
f=x;
% f=1-x*y;
%     f=1+2*x*y;
end