Tp=[1 5 10 50 100];

legendLabels="";
transformVFcFigure=figure
step(transformVFc)
hold on
for i=1:length(Tp)
    step(c2d(transformVFc,Tp(i)))
    hold on
    legendLabels(i)="Tp: "+Tp(i);
end
legend(legendLabels, 'Location','best');
title("Transmitancja objętości w zbiorniku po skoku Fc")

legendLabels="";
transformVFhFigure=figure
step(transformVFh)
hold on
for i=1:length(Tp)
    step(c2d(transformVFh,Tp(i)))
    hold on
    legendLabels(i)="Tp: "+Tp(i);
end
legend(legendLabels, 'Location','best');
title("Transmitancja objętości w zbiorniku po skoku Fh")

legendLabels="";
transformVFdFigure=figure
step(transformVFd)
hold on
for i=1:length(Tp)
    step(c2d(transformVFd,Tp(i)))
    hold on
    legendLabels(i)="Tp: "+Tp(i);
end
legend(legendLabels, 'Location','best');
title("Transmitancja objętości w zbiorniku po skoku Fd")

legendLabels="";
transformHFcFigure=figure
step(transformHFc)
hold on
for i=1:length(Tp)
    step(c2d(transformHFc,Tp(i)))
    hold on
    legendLabels(i)="Tp: "+Tp(i);
end
legend(legendLabels, 'Location','best');
title("Transmitancja wysokości słupa cieczy w zbiorniku po skoku Fc")

legendLabels="";
transformHFhFigure=figure
step(transformHFh)
hold on
for i=1:length(Tp)
    step(c2d(transformHFh,Tp(i)))
    hold on
    legendLabels(i)="Tp: "+Tp(i);
end
legend(legendLabels, 'Location','best');
title("Transmitancja wysokości słupa cieczy w zbiorniku po skoku Fh")

legendLabels="";
transformHFdFigure=figure
step(transformHFd)
hold on
for i=1:length(Tp)
    step(c2d(transformHFd,Tp(i)))
    hold on
    legendLabels(i)="Tp: "+Tp(i);
end
legend(legendLabels, 'Location','best');
title("Transmitancja wysokości słupa cieczy w zbiorniku po skoku Fd")

legendLabels="";
transformTFcFigure=figure
step(transformTFc)
hold on
for i=1:length(Tp)
    step(c2d(transformTFc,Tp(i)))
    hold on
    legendLabels(i)="Tp: "+Tp(i);
end
legend(legendLabels, 'Location','best');
title("Transmitancja temperatury w zbiorniku po skoku Fc")

legendLabels="";
transformTFhFigure=figure
step(transformTFh)
hold on
for i=1:length(Tp)
    step(c2d(transformTFh,Tp(i)))
    hold on
    legendLabels(i)="Tp: "+Tp(i);
end
legend(legendLabels, 'Location','best');
title("Transmitancja temperatury w zbiorniku po skoku Fh")

legendLabels="";
transformTFdFigure=figure
step(transformTFd)
hold on
for i=1:length(Tp)
    step(c2d(transformTFd,Tp(i)))
    hold on
    legendLabels(i)="Tp: "+Tp(i);
end
legend(legendLabels, 'Location','best');
title("Transmitancja temperatury w zbiorniku po skoku Fd")

legendLabels="";
transformTTcFigure=figure
step(transformTTc)
hold on
for i=1:length(Tp)
    step(c2d(transformTTc,Tp(i)))
    hold on
    legendLabels(i)="Tp: "+Tp(i);
end
legend(legendLabels, 'Location','best');
title("Transmitancja temperatury w zbiorniku po hipotetycznym skoku Tc")

legendLabels="";
transformTThFigure=figure
step(transformTTh)
hold on
for i=1:length(Tp)
    step(c2d(transformTTh,Tp(i)))
    hold on
    legendLabels(i)="Tp: "+Tp(i);
end
legend(legendLabels, 'Location','best');
title("Transmitancja temperatury w zbiorniku po hipotetycznym skoku Th")

legendLabels="";
transformTTdFigure=figure
step(transformTTd)
hold on
for i=1:length(Tp)
    step(c2d(transformTTd,Tp(i)))
    hold on
    legendLabels(i)="Tp: "+Tp(i);
end
legend(legendLabels, 'Location','best');
title("Transmitancja temperatury w zbiorniku po hipotetycznym skoku Td")


fileName="transforms";
overLeafFilePath="img/transforms/";
path=overLeafFilePath;


captionV="Wykresy dla transmitancji objętości";
labelV="fig:transformV";

saveFiguresInColumn([transformVFcFigure,transformVFhFigure,transformVFdFigure], path,...
    ["transformVFc","transformVFh","transformVFd"],fileName,overLeafFilePath,captionV,labelV);
    

captionH="Wykresy dla transmitancji wysokości słupa cieczy w zbiorniku";
labelH="fig:transformH";

% saveFiguresInColumn([transformHFcFigure,transformHFhFigure,transformHFdFigure], path,...
%     ["transformHFc","transformHFh","transformHFd"],fileName,overLeafFilePath,captionH,labelH);
% 
% 
% captionTF="Wykresy dla transmitancji temperatury wyjściowej po zmianie dopływu cieczy";
% labelTF="fig:transformTF";
% 
% saveFiguresInColumn([transformTFcFigure,transformTFhFigure,transformTFdFigure], path,...
%     ["transformTFc","transformTFh","transformTFd"],fileName,overLeafFilePath,captionTF,labelTF);
% 
% 
captionTT="Wykresy dla transmitancji temperatury wyjściowej po zmianie temperatury cieczy";
labelTT="fig:transformTT";

saveFiguresInColumn([transformTTcFigure,transformTThFigure,transformTTdFigure], path,...
    ["transformTTc","transformTTh","transformTTd"],fileName,overLeafFilePath,captionTT,labelTT);