% STEP RESPONSES

%init variables
init
    
Tinputs=[Th,Tc,Td];
F0inputs=[Fh,Fcin,Fd];

FcVector =Fcin +[-20,-10,0,10,20]; 
FhVector = Fh +[-14,-10,0,10,20]; 

% FcVector=Fcin+10;
FcVectorLength = length(FcVector);
FhVectorLength = length(FhVector);


hVector=ones(length(t),FcVectorLength,FhVectorLength)*h0;
hVectorLinearized=ones(length(t),FcVectorLength,FhVectorLength)*h0;

V=zeros(length(t),FcVectorLength,FhVectorLength);
V0=volume(hVector(1));
V(1,:,:)=V0;
VL=zeros(length(t),FcVectorLength,FhVectorLength);
VL(1,:,:)=volume(hVectorLinearized(1));

Tvector=ones(length(t),FcVectorLength,FhVectorLength)*T0;
TvectorL=ones(length(t),FcVectorLength,FhVectorLength)*T0;

TvectorOutput=ones(length(t),FcVectorLength,FhVectorLength)*T0;
TvectorLOutput=ones(length(t),FcVectorLength,FhVectorLength)*T0;

for j=1:length(FhVector)
    Fh=FhVector(j);
    
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
         
        h=hVector(k-1,i,j);
        hL=hVectorLinearized(k-1,i,j);
        
        T=Tvector(k-1,i,j);
        TL=TvectorL(k-1,i,j);
        
        kV1= dVdt(h, delay, Finputs);
        kV2= dVdt(heightFromVolume(V(k-1,i,j) + Tp/2*kV1),delay,Finputs);
        kV3= dVdt(heightFromVolume(V(k-1,i,j) + Tp/2*kV2),delay,Finputs);
        kV4= dVdt(heightFromVolume(V(k-1,i,j) + Tp*kV3),delay,Finputs);

        dV=Tp/6*(kV1+2*kV2+2*kV3+kV4);
        V(k,i,j)=V(k-1,i,j)+dV;
        hVector(k,i,j)=heightFromVolume(V(k,i,j));

        kV1L= dVdtLinearized(hL,h0, delay, Finputs);
        kV2L= dVdtLinearized(heightFromVolume(VL(k-1,i,j) + Tp/2*kV1L),h0,delay,Finputs);
        kV3L= dVdtLinearized(heightFromVolume(VL(k-1,i,j) + Tp/2*kV2L),h0,delay,Finputs);
        kV4L= dVdtLinearized(heightFromVolume(VL(k-1,i,j) + Tp*kV3L),h0,delay,Finputs);

        dVL=Tp/6*(kV1L+2*kV2L+2*kV3L+kV4L);
        VL(k,i,j)=VL(k-1,i,j)+dVL;
        hVectorLinearized(k,i,j)=heightFromVolume(VL(k,i,j));
        
                
        kT1= dTdt(V(k-1,i,j),T,delay,Finputs,Tinputs);
        kT2= dTdt(V(k-1,i,j) + Tp*V(k-1,i,j)/2,T + Tp/2*kT1,delay,Finputs,Tinputs);
        kT3= dTdt(V(k-1,i,j) + Tp*V(k-1,i,j)/2,T + Tp/2*kT2,delay,Finputs,Tinputs);
        kT4= dTdt(V(k-1,i,j) + Tp*V(k-1,i,j),T + Tp*kT3,delay,Finputs,Tinputs);       

        dT=Tp/6*(kT1+2*kT2+2*kT3+kT4);
        Tvector(k,i,j)=Tvector(k-1,i,j)+ dT;
        

        kT1L= dTdtLinearized(VL(k-1,i,j),V0,TL,T0,delay,Finputs,F0inputs,Tinputs);
        kT2L= dTdtLinearized(VL(k-1,i,j) + Tp/2*kV1L,V0,TL + Tp/2*kT1L,T0,delay,Finputs,F0inputs,Tinputs);
        kT3L= dTdtLinearized(VL(k-1,i,j) + Tp/2*kV2L,V0,TL + Tp/2*kT2L,T0,delay,Finputs,F0inputs,Tinputs);
        kT4L= dTdtLinearized(VL(k-1,i,j) + Tp*kV3L,V0,TL + Tp*kT3L,T0,delay,Finputs,F0inputs,Tinputs);
        
        dTLinearized=Tp/6*(kT1L+2*kT2L+2*kT3L+kT4L);
        TvectorL(k,i,j)=TvectorL(k-1,i,j)+dTLinearized;
        
        if(k<=delayT)
            TvectorOutput(k,i,j)=T0;
            TvectorLOutput(k,i,j)=T0;
        else
            TvectorOutput(k,i,j)=Tvector(k-delayT,i,j);
            TvectorLOutput(k,i,j)=TvectorL(k-delayT,i,j);
        end
    end
end
end

fileName="step-responses-h";
overLeafFilePath="img/step-responses/h/";
path=overLeafFilePath;

legendLabels=[""];
for j=1:FhVectorLength
    figure
    for i=1:FcVectorLength
        plot(t,hVector(:,i,j), 'r')
        hold on
        plot(t,hVectorLinearized(:,i,j),'g')
        hold on
        legendLabels(2*i-1)="Fc[$\frac{cm^3}{s}$]:" +FcVector(i);
        legendLabels(2*i)="linearyzacja Fc[$\frac{cm^3}{s}$]:" +FcVector(i);
    end
    title("Odpowiedzi skokowe h dla"+newline+"Fh[$\frac{cm^3}{s}$]:" +FhVector(j),'Interpreter', 'latex');
    legend(legendLabels, 'Location', 'best','Interpreter','latex')
    xlabel("t[s]");
    ylabel("h[cm]")
    hold off
    
    name="stepResponseHFh"+FhVector(j);
    caption="Poziom cieczy w zbiorniku w odpowiedzi skokowej dla skoku Fh[$\frac{cm^3}{s}$]: "+FhVector(j);
    label="fig:stepResponseHFh"+FhVector(j);

%     saveFigure(gcf, path,name,fileName,overLeafFilePath,caption,label);
end

for j=1:FhVectorLength
    figure
    for i=1:FcVectorLength
        plot(t,Tvector(:,i,j), 'r')
        hold on
        plot(t,TvectorL(:,i,j),'g')
        hold on
    end
    title("Odpowiedzi skokowe temperatura TODO "+FhVector(j))
    legend("charakterystyka dynamiczna","charakterystyka zlinearyzowana", 'Location', 'best')
    xlabel("t[s]");
    ylabel("T[\circC]")
    hold off
end

fileName="step-responses-Tout";
overLeafFilePath="img/step-responses/Tout/";
path=overLeafFilePath;

legendLabels=[""];
for j=1:FhVectorLength
    figure
    for i=1:FcVectorLength
        plot(t,TvectorOutput(:,i,j), 'r')
        hold on
        plot(t,TvectorLOutput(:,i,j),'g')
        hold on
        legendLabels(2*i-1)="Fc[$\frac{cm^3}{s}$]:" +FcVector(i);
        legendLabels(2*i)="linearyzacja Fc[$\frac{cm^3}{s}$]:" +FcVector(i);
    end
    title("Odpowiedzi skokowe temperatury wyjsciowej dla"+newline+"Fh[$\frac{cm^3}{s}$]:" +FhVector(j),'Interpreter', 'latex');
    legend(legendLabels, 'Location', 'best', 'Interpreter','latex')
    xlabel("t[s]");
    ylabel("T[\circC]")
    hold off
    
    name="stepResponseToutFh"+FhVector(j);
    caption="Temperatura na wyjściu zbiornika w odpowiedzi skokowej Fh[$"+"\\"+"frac{cm^3}{s}$]: "+FhVector(j);
    label="fig:stepResponseToutFh"+FhVector(j);

    
%     saveFigure(gcf, path,name,fileName,overLeafFilePath,caption,label);
end
