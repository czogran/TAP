% STEP RESPONSES

%init variables
init
    
Tinputs=[Th,Tc,Td];
T0inputs=[Th,Tc,Td];

F0inputs=[Fh,Fcin,Fd];

FhVector = Fh +[-14,-10,0,10,20]; 
Fc=Fcin;

FhVectorLength = length(FhVector);

FhU=ones(length(t),FhVectorLength)*Fh;
FcU=ones(length(t),1)*Fc;


hVector=ones(length(t),FhVectorLength)*h0;
hVectorLinearized=ones(length(t),FhVectorLength)*h0;

V=zeros(length(t),FhVectorLength);
V(1,:)=V0;
VL=zeros(length(t),FhVectorLength);
VL(1,:)=volume(hVectorLinearized(1));

Tvector=ones(length(t),FhVectorLength)*T0;
TvectorL=ones(length(t),FhVectorLength)*T0;

TvectorOutput=ones(length(t),FhVectorLength)*T0;
TvectorLOutput=ones(length(t),FhVectorLength)*T0;


    
    for i=1:length(FhVector)    
    Fh=FhVector(i);
    % rungy-kutta
    Finputs=[Fh,Fc,Fd];

    for k=2:length(t)
        delay=1;
        FhU(k,i)=Fh;
         
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

fileName="step-responses-fh";
overLeafFilePath="img/step-responses-fh/";
path=overLeafFilePath;

    heightFigure=figure
    for i=1:FhVectorLength
        plot(t,hVector(:,i),'Color',colorLabels(i))
        hold on
        plot(t,hVectorLinearized(:,i),'Color',colorLabels(i),'LineStyle',"--")
        hold on
    end
    title("Odpowiedzi skokowe h dla Fh przy"+newline+"Fc[$\frac{cm^3}{s}$]:" +Fcin,'Interpreter', 'latex');
    xlabel("t[s]");
    ylabel("h[cm]")
    hold off
    
    heightName="stepResponseHFh";
     
    legendLabels=[""];
    tempFigure=figure
    for i=1:FhVectorLength
        plot(t,TvectorOutput(:,i),'Color',colorLabels(i))
        hold on
        plot(t,TvectorLOutput(:,i),'Color',colorLabels(i),'LineStyle',"--")
        hold on
    end
    title("Odpowiedzi skokowe temperatury wyjsciowej dla Fh przy"+newline+"Fc[$\frac{cm^3}{s}$]:" +Fcin,'Interpreter', 'latex');
    xlabel("t[s]");
    ylabel("T[\circC]")
    hold off
    
    tempName="stepResponseToutFh";

    legendLabels="";
    controlFigure=figure
    sgtitle('Sterowania Fh i Fc') 
    
    subplot(2,1,1)
    plot(t,FcU,'Color',colorLabels(i));
    xlabel("t[s]");
    ylabel("Fc[$\frac{cm^3}{s}$]",'Interpreter', 'latex')
    title("Sterowanie Fc")
   
    axis([-50 t(end) min(FcU(:,1))-5 max(FcU(:,1))+5])
    hold off
    
    subplot(2,1,2)
    for i=1:FhVectorLength
        plot(t,FhU(:,i),'Color',colorLabels(i))
        hold on
    end
    title("Sterowanie Fh");
    xlabel("t[s]");
    ylabel("Fh[$\frac{cm^3}{s}$]",'Interpreter', 'latex')
    hold off
    
    controlName="stepResponseUFh";
    
    caption="Wykresy dla odpowiedzi skokowej Fh przy Fc[$"+"\\"+"frac{cm^3}{s}$]: "+Fcin;
    label="fig:stepResponsesFh";

    
    legendName="legendFc";

legendLabels="";
legendFigure=figure
for i=1:FhVectorLength
    plot(0,'Color',colorLabels(i))
    hold on
    plot(0,'Color',colorLabels(i),'LineStyle',"--")
    hold on
    legendLabels(2*i-1)="Fh[$\frac{cm^3}{s}$]:" +FhVector(i);
    legendLabels(2*i)="linearyzacja Fh[$\frac{cm^3}{s}$]:" +FhVector(i);
end
title("Legenda",'Interpreter', 'latex');
legend(legendLabels, 'Location', 'best', 'Interpreter','latex', 'FontSize',12)
hold off
    saveFiguresInColumn([legendFigure,heightFigure,tempFigure,controlFigure], path,[legendName,heightName,tempName,controlName],fileName,overLeafFilePath,caption,label);

