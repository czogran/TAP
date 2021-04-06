transformVFcFigure=figure
step(transformVFc)
title("Transmitancja objętości w zbiorniku po skoku Fc")

transformVFhFigure=figure
step(transformVFh)
title("Transmitancja objętości w zbiorniku po skoku Fh")

transformVFdFigure=figure
step(transformVFd)
title("Transmitancja objętości w zbiorniku po skoku Fd")


transformHFcFigure=figure
step(transformHFc)
title("Transmitancja wysokości słupa cieczy w zbiorniku po skoku Fc")

transformHFhFigure=figure
step(transformHFh)
title("Transmitancja wysokości słupa cieczy w zbiorniku po skoku Fh")

transformHFdFigure=figure
step(transformHFd)
title("Transmitancja wysokości słupa cieczy w zbiorniku po skoku Fd")

transformTFcFigure=figure
step(transformTFc)
title("Transmitancja temperatury w zbiorniku po skoku Fc")

transformTFhFigure=figure
step(transformTFh)
title("Transmitancja temperatury w zbiorniku po skoku Fh")

transformTFdFigure=figure
step(transformTFd)
title("Transmitancja temperatury w zbiorniku po skoku Fd")


transformTTcFigure=figure
step(transformTTc)
title("Transmitancja temperatury w zbiorniku po hipotetycznym skoku Tc")

transformTThFigure=figure
step(transformTTh)
title("Transmitancja temperatury w zbiorniku po hipotetycznym skoku Th")

transformTFTFigure=figure
step(transformTTd)
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

saveFiguresInColumn([transformHFcFigure,transformHFhFigure,transformHFdFigure], path,...
    ["transformHFc","transformHFh","transformHFd"],fileName,overLeafFilePath,captionH,labelH);


captionTF="Wykresy dla transmitancji temperatury wyjściowej po zmianie dopływu cieczy";
labelTF="fig:transformTF";

saveFiguresInColumn([transformTFcFigure,transformTFhFigure,transformTFdFigure], path,...
    ["transformTFc","transformTFh","transformTFd"],fileName,overLeafFilePath,captionTF,labelTF);


captionTT="Wykresy dla transmitancji temperatury wyjściowej po zmianie temperatury cieczy";
labelTT="fig:transformTT";

saveFiguresInColumn([transformTTcFigure,transformTThFigure,transformTThFigure], path,...
    ["transformTTc","transformTTh","transformTTh"],fileName,overLeafFilePath,captionTT,labelTT);