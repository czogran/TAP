a=0.25;
a1=0.2;

N=3;
Nu=2;

yPrev=[2 2 1];
uPrev=[0.5 1];
uk=[1.35 1.77 1.2668];
u=[uPrev,uk ];

y=[yPrev,zeros(1,N)];

d=yPrev(end)-(a*yPrev(end-1)^2+uPrev(end)+a1*yPrev(end-2)*uPrev(end-1));

for k=2:2+N
    y(k+1)=a*y(k)^2+u(k)+a1*y(k-1)*u(k-1)+d;
end

y