% STEP RESPONSES

%init variables
init
    
Tinputs=[Th,Tc,Td];
T0inputs=[Th,Tc,Td];

F0inputs=[Fh,Fcin,Fd];

FcVector =Fcin +[-20,-10,0,10,20]; 

FcVectorLength = length(FcVector);

FcU=ones(length(t),FcVectorLength)*Fcin;
FhU=ones(length(t),FcVectorLength)*Fh;

hVector=ones(length(t),FcVectorLength)*h0;
hVectorLinearized=ones(length(t),FcVectorLength)*h0;

V=zeros(length(t),FcVectorLength);
V(1,:,:)=V0;
VL=zeros(length(t),FcVectorLength);
VL(1,:,:)=volume(hVectorLinearized(1));

Tvector=ones(length(t),FcVectorLength)*T0;
TvectorL=ones(length(t),FcVectorLength)*T0;

TvectorOutput=ones(length(t),FcVectorLength)*T0;
TvectorLOutput=ones(length(t),FcVectorLength)*T0;

for i=1:length(FcVector)    
% rungy-kutta
for k=2:length(t)
    delay=1;

    if(k<delayC)
        Fc=Fcin;
    else
       Fc=FcVector(i);
    end 
    Finputs=[Fh,Fc,Fd];
    FhU(k,i)=Fh;
    FcU(k,i)=Fc;


    h=hVector(k-1,i);
    hL=hVectorLinearized(k-1,i);

    T=Tvector(k-1,i);
    TL=TvectorL(k-1,i);

    kV1= dVdt(h, delay, Finputs);
    kV2= dVdt(heightFromVolume(V(k-1,i) + Tp/2*kV1),delay,Finputs);
    kV3= dVdt(heightFromVolume(V(k-1,i) + Tp/2*kV2),delay,Finputs);
    kV4= dVdt(heightFromVolume(V(k-1,i) + Tp*kV3),delay,Finputs);

    dV=Tp/6*(kV1+2*kV2+2*kV3+kV4);
    V(k,i)=V(k-1,i)+dV;
    hVector(k,i)=heightFromVolume(V(k,i));

    kV1L= dVdtLinearized(hL,h0, delay, Finputs);
    kV2L= dVdtLinearized(heightFromVolumeLinearized(VL(k-1,i) + Tp/2*kV1L,V0),h0,delay,Finputs);
    kV3L= dVdtLinearized(heightFromVolumeLinearized(VL(k-1,i) + Tp/2*kV2L,V0),h0,delay,Finputs);
    kV4L= dVdtLinearized(heightFromVolumeLinearized(VL(k-1,i) + Tp*kV3L,V0),h0,delay,Finputs);

    dVL=Tp/6*(kV1L+2*kV2L+2*kV3L+kV4L);
    VL(k,i)=VL(k-1,i)+dVL;
    hVectorLinearized(k,i)=heightFromVolumeLinearized(VL(k,i),V0);


    kT1= dTdt(V(k-1,i),T,delay,Finputs,Tinputs);
    kT2= dTdt(V(k-1,i) + Tp*V(k-1,i)/2,T + Tp/2*kT1,delay,Finputs,Tinputs);
    kT3= dTdt(V(k-1,i) + Tp*V(k-1,i)/2,T + Tp/2*kT2,delay,Finputs,Tinputs);
    kT4= dTdt(V(k-1,i) + Tp*V(k-1,i),T + Tp*kT3,delay,Finputs,Tinputs);       

    dT=Tp/6*(kT1+2*kT2+2*kT3+kT4);
    Tvector(k,i)=Tvector(k-1,i)+ dT;


    kT1L= dTdtLinearized(VL(k-1,i),V0,TL,T0,delay,Finputs,F0inputs,Tinputs,T0inputs);
    kT2L= dTdtLinearized(VL(k-1,i) + Tp/2*kV1L,V0,TL + Tp/2*kT1L,T0,delay,Finputs,F0inputs,Tinputs,T0inputs);
    kT3L= dTdtLinearized(VL(k-1,i) + Tp/2*kV2L,V0,TL + Tp/2*kT2L,T0,delay,Finputs,F0inputs,Tinputs,T0inputs);
    kT4L= dTdtLinearized(VL(k-1,i) + Tp*kV3L,V0,TL + Tp*kT3L,T0,delay,Finputs,F0inputs,Tinputs,T0inputs);

    dTLinearized=Tp/6*(kT1L+2*kT2L+2*kT3L+kT4L);
    TvectorL(k,i)=TvectorL(k-1,i)+dTLinearized;

    if(k<=delayT)
        TvectorOutput(k,i)=T0;
        TvectorLOutput(k,i)=T0;
    else
        TvectorOutput(k,i)=Tvector(k-delayT,i);
        TvectorLOutput(k,i)=TvectorL(k-delayT,i);
    end
end
end

fileName="step-responses-fc";
overLeafFilePath="img/step-responses-fc/";
path=overLeafFilePath;

heightFigure=figure
for i=1:FcVectorLength
    plot(t,hVector(:,i),'Color',colorLabels(i))
    hold on
    plot(t,hVectorLinearized(:,i),'Color',colorLabels(i),'LineStyle',"--")
    hold on
end
title("Odpowiedzi skokowe h dla Fc przy"+newline+"Fh[$\frac{cm^3}{s}$]:" +Fh,'Interpreter', 'latex');
xlabel("t[s]");
ylabel("h[cm]")
hold off

heightName="stepResponseHFc"

tempFigure=figure
for i=1:FcVectorLength
    plot(t,TvectorOutput(:,i),'Color',colorLabels(i))
    hold on
    plot(t,TvectorLOutput(:,i),'Color',colorLabels(i),'LineStyle',"--")
    hold on
end
title("Odpowiedzi skokowe temperatury wyjsciowej dla Fc przy"+newline+"Fh[$\frac{cm^3}{s}$]:" +Fh,'Interpreter', 'latex');
xlabel("t[s]");
ylabel("T[\circC]")
hold off

tempName="stepResponseToutFh";

legendLabels="";
controlFigure=figure
sgtitle('Sterowania Fh i Fc') 

subplot(2,1,1)
plot(t,FhU(:,1),'Color',colorLabels(i));
xlabel("t[s]");
ylabel("Fh[$\frac{cm^3}{s}$]",'Interpreter', 'latex')
title("Sterowanie Fh")

axis([-50 t(end) min(FhU(:,1))-5 max(FhU(:,1))+5])
hold off

subplot(2,1,2)
for i=1:FcVectorLength
    plot(t,FcU(:,i),'Color',colorLabels(i))
    hold on
    legendLabels(i)="Fc sterowanie: " +FcVector(i);
end
title("Sterowanie Fc");
xlabel("t[s]");
ylabel("Fc[$\frac{cm^3}{s}$]",'Interpreter', 'latex')
hold off

controlName="stepResponseUFc";




legendName="legendFc";

legendLabels="";
legendFigure=figure
for i=1:FcVectorLength
    plot(0,'Color',colorLabels(i))
    hold on
    plot(0,'Color',colorLabels(i),'LineStyle',"--")
    hold on
    legendLabels(2*i-1)="Fc[$\frac{cm^3}{s}$]:" +FcVector(i);
    legendLabels(2*i)="linearyzacja Fc[$\frac{cm^3}{s}$]:" +FcVector(i);
end
title("Legenda",'Interpreter', 'latex');
legend(legendLabels, 'Location', 'best', 'Interpreter','latex', 'FontSize',12)
hold off


caption="Wykresy dla odpowiedzi skokowej Fc dla Fh[$"+"\\"+"frac{cm^3}{s}$]: "+Fh;
label="fig:stepResponsesFc"+Fh;
saveFiguresInColumn([legendFigure,heightFigure,tempFigure,controlFigure], path,[legendName,heightName,tempName,controlName],fileName,overLeafFilePath,caption,label);

