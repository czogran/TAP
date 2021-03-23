% EXAMPLE USAGE
% dodanie �cie�ki
    % w skrypcie nale�y doda� �cie�k� addpath('..\common\');
% deklaracja figury
    % figure
% jakie� ploty
    % plot(y) 
% gdy warto�ci maj� kropki trzeba je zamieni�, by m�c zapisa� w pliku
    % Estr=sprintf('%0.5f',E);
    % Estr=strrep(Estr,'.','-' );
%  nazwa wykresu wa�ne by pisa� w "" a nie w '' bo wtdy nie dzia�a, robienie stringa w
%  zale�no�ci od parametr�w wykresu
    % nazwa = "DMC_d"+D+"_N"+N+"_Nu"+Nu+"_lam"+lambda+"_e"+Estr; 
% podpis jaki b�dzie pod wykresem na overleafie
    % caption="Wykres $D=\\num{"+D+"},\\  N = \\num{"+N+"},\\  N_{\\mathrm{u}} = \\num{"+Nu+"}, \\lambda =\\num{"+lambda+"},\\  E=\\num{"+round(E,3)+"}$";
% labelka w overleafie, warto da� ja od parametru kt�ry si� zmienia, wtedy �atwo si� mo�na do niej dosta� 
    % label="fig:lambda"+lambda;
%  nazwa dokumentu w kt�rym b�dzie skrypt do rysowania obrazk�w
    % fileName="lambda";
% �cie�ka w kt�rej b�d� si� na overleafie rysunki znajdowa�
    % overLeafFilePath="img/5_img/DMC/lambda/";
% zapis do pliku     
    % saveFigure(gcf, 'img/5img/DMC/',nazwa,fileName,overLeafFilePath,caption,label);

% WA�NA UWAGA, by nadpisa� obecne rzeczy w pliku nale�y odpali� skrypt po
% 10 sod ostatniego wywo�ania, jest timer kt�ry kasuje wtedy zawartos�
% pliku
function  saveFigure(figure,path,name,fileName,overLeafFilePath,caption,label)
    saveas(figure,path+name,'epsc')
    
 
    
    fileContent={
    '\\begin{figure}[tb]'
    '   \\centering'
    '   \\includegraphics{'+overLeafFilePath+name+'.eps}'
    '   \\caption{'+caption+'}'
    '   \\label{'+label+'}'
    '\\end{figure}'
    '            '
    };
    fileContent=sprintf('%s\n',fileContent{:})
    
    global previousUsageTime;
    if (isempty(previousUsageTime))
        fileID = fopen(path+fileName+'.tex','w')
    else
        disp(posixtime(datetime())-previousUsageTime);
        if ((posixtime(datetime())-previousUsageTime)>10)
            fileID = fopen(path+fileName+'.tex','w');
        else
            
             fileID = fopen(path+fileName+'.tex','a');
        end
    end
    previousUsageTime=posixtime(datetime());

    fprintf(fileID,fileContent);
    fclose(fileID);
    
end

