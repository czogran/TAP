% STEP RESPONSES

%init variables
init

FcVector =Fcin +[-20,-10,0,10,20]; 
% FcVector=Fcin+10;
FcVectorLength = length(FcVector);

hVector=ones(length(t),FcVectorLength)*h0;
hVectorLinearized=ones(length(t),FcVectorLength)*h0;

V=zeros(length(t),FcVectorLength);
V(1,:)=volume(hVector(1));
VL=zeros(length(t),FcVectorLength);
VL(1,:)=volume(hVectorLinearized(1));

Tvector=[0];

for i=1:length(FcVector)    
    Fc=FcVector(i);
    Finputs=[Fh,Fc,Fd];

    % rungy-kutta
    for k=2:length(t)
        h=hVector(k-1,i);
        hL=hVectorLinearized(k-1,i);

        delay=1;

        kV1= dVdt(h, delay, Finputs);
        kV2= dVdt(heightFromVolume(V(k-1,i) + Tp/2*kV1),delay,Finputs);
        kV3= dVdt(heightFromVolume(V(k-1,i) + Tp/2*kV2),delay,Finputs);
        kV4= dVdt(heightFromVolume(V(k-1,i) + Tp*kV3),delay,Finputs);

        dV=Tp/6*(kV1+2*kV2+2*kV3+kV4);
        V(k,i)=V(k-1,i)+dV;
        hVector(k,i)=heightFromVolume(V(k,i));

        kV1L= dVdtLinearized(hL,h0, delay, Finputs);
        kV2L= dVdtLinearized(heightFromVolume(VL(k-1,i) + Tp/2*kV1L),h0,delay,Finputs);
        kV3L= dVdtLinearized(heightFromVolume(VL(k-1,i) + Tp/2*kV2L),h0,delay,Finputs);
        kV4L= dVdtLinearized(heightFromVolume(VL(k-1,i) + Tp*kV3L),h0,delay,Finputs);

        dVL=Tp/6*(kV1L+2*kV2L+2*kV3L+kV4L);
        VL(k,i)=VL(k-1,i)+dVL;
        hVectorLinearized(k,i)=heightFromVolume(VL(k,i));
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