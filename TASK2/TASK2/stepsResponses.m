% STEP RESPONSES

%init variables
init

FcVector =Fcin +[-20,-10,0,10,20]; 
FcVectorLength = length(FcVector);

hVector=ones(length(t),FcVectorLength)*h0;
hVectorLinearized=ones(length(t),FcVectorLength)*h0;

V=zeros(length(t),FcVectorLength);
VLinearized=zeros(length(t),FcVectorLength);
V(1,:)=volume(hVector(1));
VLinearized(1,:)=volume(hVectorLinearized(1));

Tvector=[0];

for i=1:length(FcVector)
    
    Fc=FcVector(i);

    Finputs=[Fh,Fc,Fd];

    % rungy-kutta
    for k=2:length(t)
        h=hVector(k-1,i);
        hL=hVectorLinearized(k-1,i);

    %     if(k<delayC)
    %       delay=0;  
    %     else
    %        delay=1;
    %     end

        delay=1;

        k1= dVdt(h, delay, Finputs);
        k2= dVdt(h+Tp/2*k1,delay,Finputs);
        k3= dVdt(h + Tp/2*k2,delay,Finputs);
        k4= dVdt(h+Tp*k3,delay,Finputs);

        dV=Tp/6*(k1+2*k2+2*k3+k4);
        V(k,i)=V(k-1,i)+dV;
        hVector(k,i)=heightFromVolume(V(k,i));

        k1L= dVdtLinearized(hL,h0, delay, Finputs);
        k2L= dVdtLinearized(hL+Tp/2*k1,h0,delay,Finputs);
        k3L= dVdtLinearized(hL + Tp/2*k2,h0,delay,Finputs);
        k4L= dVdtLinearized(hL+Tp*k3,h0,delay,Finputs);

        dVLinearized=Tp/6*(k1L+2*k2L+2*k3L+k4L);
    %     dVLinearized=dVdtLinearized(hL,h0, delay, Finputs);

        VLinearized(k,i)=VLinearized(k-1,i)+dVLinearized;
        hVectorLinearized(k,i)=heightFromVolume(VLinearized(k,i));

    %     Tvector(k)=Tvector(k-1)+ dVdTdt/V(k);
    end
end

figure
for i=1:FcVectorLength
    plot(t,hVector(:,i), 'r')
    hold on
    plot(t,hVectorLinearized(:,i),'g')
    hold on
end
title("Odpowiedzi skokowe TODO")
legend("charakterystyka dynamiczna","charakterystyka zlinearyzowana", 'Location', 'best')
xlabel("t[s]");
ylabel("h[cm]")
hold off

%     dVLinearized=Tp*dVdtLinearized(hLinearized,h0,delay, Finputs);
%     VLinearized(k)=V(k-1)+dV;
%     hVector(k)=heightFromVolume(V(k));
% %     Tvector(k)=Tvector(k-1)+ dVdTdt/V(k);