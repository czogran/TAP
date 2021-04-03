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


Tvector=[0];
Tinputs=[Th,Tc,Td];

Fc=Fcin;

FhStatic = (Fh-10) : (Fh+10);
FcStatic = (Fc-10) : (Fc+10);

Vstatic = zeros(length(FhStatic), length(FcStatic));
hStatic = zeros(length(FhStatic), length(FcStatic));
Tstatic = zeros(length(FhStatic), length(FcStatic));

for i =  1: length(FhStatic)
    Fh=FhStatic(i);
    for  j=  1: length(FcStatic)
        Fc=FcStatic(j);
        h=0;
        V=[0];
        Tvector=[0];
        
        Finputs=[Fh, Fc,Fd];
        for k=2:length(t)
            h=heightFromVolume(V(k-1));
            T=Tvector(k-1);
            
            kV1= dVdt(h, delay, Finputs);
            kV2= dVdt(heightFromVolume(V(k-1) + Tp/2*kV1),delay,Finputs);
            kV3= dVdt(heightFromVolume(V(k-1) + Tp/2*kV2),delay,Finputs);
            kV4= dVdt(heightFromVolume(V(k-1) + Tp*kV3),delay,Finputs);


            dV=Tp/6*(kV1+2*kV2+2*kV3+kV4);
            V(k)=V(k-1)+dV;
            hVector(k)=heightFromVolume(V(k));
            
            
            kT1= dTdt(V(k-1),T,delay,Finputs,Tinputs);
            kT2= dTdt(V(k-1) + Tp*V(k-1)/2,T + Tp/2*kT1,delay,Finputs,Tinputs);
            kT3= dTdt(V(k-1) + Tp*V(k-1)/2,T + Tp/2*kT2,delay,Finputs,Tinputs);
            kT4= dTdt(V(k-1) + Tp*V(k-1),T + Tp*kT3,delay,Finputs,Tinputs);       

            dT=Tp/6*(kT1+2*kT2+2*kT3+kT4);
            Tvector(k)=Tvector(k-1)+ dT;

        end
        Vstatic(i,j) = V(end);
        hStatic(i,j)=heightFromVolume(V(end));
        Tstatic(i,j)=Tvector(end);
    end
end

VFigure=figure
surf(FcStatic, FhStatic,  hStatic)
xlabel('Fc[$\frac{cm^3}{s}$]','Interpreter', 'latex');
ylabel('Fh[$\frac{cm^3}{s}$]','Interpreter', 'latex');
zlabel('h[cm]')
title("Charakterystyka statyczna"+ newline+"wysokości słupa cieczy w zależności od"+newline+"dopływu cieplej i zimnej wody");
hold off

hFigure=figure
surf(FcStatic, FhStatic,  Vstatic)
xlabel('Fc[$\frac{cm^3}{s}$]','Interpreter', 'latex');
ylabel('Fh[$\frac{cm^3}{s}$]','Interpreter', 'latex');
zlabel('V[[${cm^3}$]','Interpreter', 'latex')
title("Charakterystyka statyczna"+ newline+"objętość cieczy w zależności od"+newline+"dopływu cieplej i zimnej wody");
hold off

TFigure=figure
surf(FcStatic, FhStatic,  Tstatic)
xlabel('Fc[$\frac{cm^3}{s}$]','Interpreter', 'latex');
ylabel('Fh[$\frac{cm^3}{s}$]','Interpreter', 'latex');
zlabel('T[\circC]')
title("Charakterystyka statyczna"+ newline+"temperatura w zbiorniku lub wyjściowa" +newline+ "w zależności od"+newline+"dopływu cieplej i zimnej wody");
hold off


fileName="static-characteristics";
overLeafFilePath="img/static/";
path=overLeafFilePath;

name="staticV";
caption="Charakterystyka statyczna pojemności zbiornika";
label="fig:staticV";

saveFigure(VFigure, path,name,fileName,overLeafFilePath,caption,label);

name="staticH";
caption="Charakterystyka statyczna wysokości słupa cieczy w zbiorniku";
label="fig:staticH";

saveFigure(hFigure, path,name,fileName,overLeafFilePath,caption,label);

name="staticT";
caption="Charakterystyka statyczna temperatury cieczy w zbiorniku, czy na wyjściu zbiornika";
label="fig:staticT";

saveFigure(TFigure, path,name,fileName,overLeafFilePath,caption,label);