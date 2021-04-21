  
init

% consts
a=7;
C0=0.7;

% sampe time
Tp=1; 

addconst = 0;

%transmit = tf({1 1 1}, {[1, a/(4*sqrt(sqrt(h0^3*C0)))] [1, a/(4*sqrt(sqrt(h0^3*C0)))] [1, a/(4*sqrt(sqrt(h0^3*C0)))]}, 'InputDelay', [0; 120; 0]);
h0 = 81;

if addconst == 1
    B = [1, 1, 1, 0, 0, 0; Th0/V0 - T0/V0, Tc0/V0 - T0/V0, Td0/V0 - T0/V0, Fh0/V0, Fc0/V0, Fd0/V0;0, 0, 0, 0, 0, 0];
    A = [-a*0.25/(sqrt(sqrt((V0^3)*C0))), 0, 1; Fh0*T0/(V0^2) + Fc0*T0/(V0^2) + Fd0*T0/(V0^2) - Fh0*Th0/(V0^2) - Fc0*Tc0/(V0^2) - Fd0*Td0/(V0^2), -(Fh0/V0 + Fc0/V0 + Fd0/V0), 0; 0, 0, 0];
    C = [1, 0, 0;0, 1, 0; 0, 0, 1];
    D = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0 ,0 , 0; 0, 0, 0, 0, 0, 0];
else
    B = [1, 1, 1, 0, 0, 0; Th0/V0 - T0/V0, Tc0/V0 - T0/V0, Td0/V0 - T0/V0, Fh0/V0, Fc0/V0, Fd0/V0];
    A = [-a*0.25/(sqrt(sqrt((V0^3)*C0))), 0; Fh0*T0/(V0^2) + Fc0*T0/(V0^2) + Fd0*T0/(V0^2) - Fh0*Th0/(V0^2) - Fc0*Tc0/(V0^2) - Fd0*Td0/(V0^2), -(Fh0/V0 + Fc0/V0 + Fd0/V0)];
    C = [1, 0;0, 1];
    D = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0 ,0 , 0];
end




transmit = ss(A, B, C, D, 'InputDelay', [0, 180, 0, 0, 180, 0], 'OutputDelay',[0, 120]);
transmit.InputName={'Fh', 'Fc','Fd', 'Th','Tc','Td'};
transmit.OutputName={'V', 'T'};


Tp=[1 5 10 50 100];
% discreteTransmit = [];
stepsFhV = {};
stepsFcV = {};
stepsFdV = {};

stepsFhT = {};
stepsFcT = {};
stepsFdT = {};

stepsThT = {};
stepsTcT = {};
stepsTdT = {};

figure

for i=1:length(Tp)
    discreteTransmit = c2d(transmit,Tp(i));
    
       step(discreteTransmit,1:Tp(i):max(t));
     hold on
     
    steps = step(discreteTransmit,1:Tp(i):max(t));
    stepsFhV{end+1}=steps(:,1,1);
    stepsFcV{end+1}=steps(:,1,2);
    stepsFdV{end+1}=steps(:,1,3);

    stepsFhT{end+1}=steps(:,2,1);
    stepsFcT{end+1}=steps(:,2,2);
    stepsFdT{end+1}=steps(:,2,3);
    
    stepsThT{end+1}=steps(:,2,4);
    stepsTcT{end+1}=steps(:,2,5);
    stepsTdT{end+1}=steps(:,2,6);
end
hold off



legendLabels="";
FhV=figure
for i=1:length(Tp)
    stairs(1:Tp(i):max(t),stepsFhV{i});
    legendLabels(i)="Tp: "+Tp(i);
    hold on
end
ylabel('V[cm^3]');
xlabel('t[s]','Interpreter', 'latex');legend(legendLabels,'Location','best');
title("Model dyskretny" + newline+ "Symlucja skoku jednostkowego F_h na objętość");
hold off 


legendLabels="";
FcV=figure
for i=1:length(Tp)
    stairs(1:Tp(i):max(t),stepsFcV{i});
    legendLabels(i)="Tp: "+Tp(i);
    hold on
end
ylabel('V[cm^3]');
xlabel('t[s]','Interpreter', 'latex');legend(legendLabels,'Location','best');
title("Model dyskretny" + newline+ "Symlucja skoku jednostkowego F_c na objętość");
hold off 

legendLabels="";
FdV=figure
for i=1:length(Tp)
    stairs(1:Tp(i):max(t),stepsFdV{i});
    legendLabels(i)="Tp: "+Tp(i);
    hold on
end
ylabel('V[cm^3]');
xlabel('t[s]','Interpreter', 'latex');legend(legendLabels,'Location','best');
title("Model dyskretny" + newline+ "Symlucja skoku jednostkowego F_d na objętość");
hold off 


legendLabels="";
FhH=figure
for i=1:length(Tp)
    s=stepsFhV{i};
    h=[]
    for k=1:length(s)
        h(k)=heightFromVolumeLinearized(V0+s(k),V0)-h0;
    end
    stairs(1:Tp(i):max(t),h);
    legendLabels(i)="Tp: "+Tp(i);
    hold on
end
ylabel('h[cm]','Interpreter', 'latex');
xlabel('t[s]','Interpreter', 'latex');legend(legendLabels,'Location','best');
title("Model dyskretny" + newline+ "Symlucja skoku jednostkowego F_h na wysokość");
hold off 


legendLabels="";
FcH=figure
for i=1:length(Tp)
    s=stepsFcV{i};
    h=[]
    for k=1:length(s)
        h(k)=heightFromVolumeLinearized(V0+s(k),V0)-h0;
    end
    stairs(1:Tp(i):max(t),h);
    legendLabels(i)="Tp: "+Tp(i);
    hold on
end
ylabel('h[cm]','Interpreter', 'latex');
xlabel('t[s]','Interpreter', 'latex');legend(legendLabels,'Location','best');
title("Model dyskretny" + newline+ "Symlucja skoku jednostkowego F_c na wysokość");
hold off 


legendLabels="";
FdH=figure
for i=1:length(Tp)
    s=stepsFdV{i};
    h=[]
    for k=1:length(s)
        h(k)=heightFromVolumeLinearized(V0+s(k),V0)-h0;
    end
    stairs(1:Tp(i):max(t),h);
    legendLabels(i)="Tp: "+Tp(i);
    hold on
end
hold off
ylabel('h[cm]','Interpreter', 'latex');
xlabel('t[s]','Interpreter', 'latex');legend(legendLabels,'Location','best');
title("Model dyskretny" + newline+ "Symlucja skoku jednostkowego F_d na wysokość");
hold off 




legendLabels="";
FhT=figure
for i=1:length(Tp)
    stairs(1:Tp(i):max(t),stepsFhT{i});
    legendLabels(i)="Tp: "+Tp(i);
    hold on
end
ylabel('T[\circC]');
xlabel('t[s]','Interpreter', 'latex');legend(legendLabels,'Location','best');
title("Model dyskretny" + newline+ "Symlucja skoku jednostkowego F_h na temperaturę wyjściową");
hold off 



legendLabels="";
FcT=figure
for i=1:length(Tp)
    stairs(1:Tp(i):max(t),stepsFcT{i});
    legendLabels(i)="Tp: "+Tp(i);
    hold on
end
ylabel('T[\circC]');
xlabel('t[s]','Interpreter', 'latex');legend(legendLabels,'Location','best');
title("Model dyskretny" + newline+ "Symlucja skoku jednostkowego F_c na temperaturę wyjściową");
hold off 


legendLabels="";
FdT=figure
for i=1:length(Tp)
    stairs(1:Tp(i):max(t),stepsFdT{i});
    legendLabels(i)="Tp: "+Tp(i);
    hold on
end
ylabel('T[\circC]');
xlabel('t[s]','Interpreter', 'latex');legend(legendLabels,'Location','best');
title("Model dyskretny" + newline+ "Symlucja skoku jednostkowego F_d na temperaturę wyjściową");
hold off 




legendLabels="";
ThT=figure
for i=1:length(Tp)
    stairs(1:Tp(i):max(t),stepsThT{i});
    legendLabels(i)="Tp: "+Tp(i);
    hold on
end
ylabel('T[\circC]');
xlabel('t[s]','Interpreter', 'latex');legend(legendLabels,'Location','best');
title("Model dyskretny" + newline+ "Symlucja skoku jednostkowego T_h na temperaturę wyjściową");
hold off 


legendLabels="";
TcT=figure
for i=1:length(Tp)
    stairs(1:Tp(i):max(t),stepsTcT{i});
    legendLabels(i)="Tp: "+Tp(i);
    hold on
end
ylabel('T[\circC]');
xlabel('t[s]','Interpreter', 'latex');legend(legendLabels,'Location','best');
title("Model dyskretny" + newline+ "Symlucja skoku jednostkowego T_c na temperaturę wyjściową");
hold off 


legendLabels="";
TdT=figure
for i=1:length(Tp)
    stairs(1:Tp(i):max(t),stepsTdT{i});
    legendLabels(i)="Tp: "+Tp(i);
    hold on
end
ylabel('T[\circC]');
xlabel('t[s]','Interpreter', 'latex');legend(legendLabels,'Location','best');
title("Model dyskretny" + newline+ "Symlucja skoku jednostkowego T_d na temperaturę wyjściową");
hold off 


fileName="discrete-step-responses";
overLeafFilePath="img/discrete-step-responses/";
path=overLeafFilePath;

saveFigure(FhV, path,"FhV",fileName,overLeafFilePath,"","");
saveFigure(FcV, path,"FcV",fileName,overLeafFilePath,"","");
saveFigure(FdV, path,"FdV",fileName,overLeafFilePath,"","");

saveFigure(FhH, path,"FhH",fileName,overLeafFilePath,"","");
saveFigure(FcH, path,"FcH",fileName,overLeafFilePath,"","");
saveFigure(FdH, path,"FdH",fileName,overLeafFilePath,"","");


saveFigure(FhT, path,"FhT",fileName,overLeafFilePath,"","");
saveFigure(FcT, path,"FcT",fileName,overLeafFilePath,"","");
saveFigure(FdT, path,"FdT",fileName,overLeafFilePath,"","");

saveFigure(ThT, path,"ThT",fileName,overLeafFilePath,"","");
saveFigure(TcT, path,"TcT",fileName,overLeafFilePath,"","");
saveFigure(TdT, path,"TdT",fileName,overLeafFilePath,"","");