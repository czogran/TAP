%% STATIC CHARACTERISTICS %%

%init data
init

%VARIABLES
% overrides init
% sample time
Tp=1;
t=0:Tp:2000;

%% static characteristics
%!!!!! IN STATIC CHARACTERISTIC DELAY DOES NOT MATTER
delay=1;

Tvector=[T];

Fc=Fcin;

FhStatic = (Fh-10) : (Fh+10);
FcStatic = (Fc-10) : (Fc+10);

Vstatic = zeros(length(FhStatic), length(FcStatic));
hStatic = zeros(length(FhStatic), length(FcStatic));


for i =  1: length(FhStatic)
    Fh=FhStatic(i);
    for  j=  1: length(FcStatic)
        Fc=FcStatic(j);
        h=81;
        V=[volume(h)];
        Tvector=[T];
        Finputs=[Fh, Fc,Fd];
        for k=2:length(t)
            kV1= dVdt(h, delay, Finputs);
            kV2= dVdt(heightFromVolume(V(k-1) + Tp/2*kV1),delay,Finputs);
            kV3= dVdt(heightFromVolume(V(k-1) + Tp/2*kV2),delay,Finputs);
            kV4= dVdt(heightFromVolume(V(k-1) + Tp*kV3),delay,Finputs);


            dV=Tp/6*(kV1+2*kV2+2*kV3+kV4);
            V(k)=V(k-1)+dV;
            h=heightFromVolume(V(k));

        end
        Vstatic(i,j) = V(end);
        hStatic(i,j)=heightFromVolume(V(end));
    end
end

figure
surf(FcStatic, FhStatic,  hStatic)
xlabel('Fh');
ylabel('Fc');
title("Charakterystyka statyczna"+ newline+"wysokości słupa cieczy w zależności od"+newline+"dopływu cieplej i zimnej wody");
hold off