% EXAMPLE USAGE
% dodanie œcie¿ki
    % w skrypcie nale¿y dodaæ œcie¿kê addpath('..\common\');
% deklaracja figury
    % figure
% jakieœ ploty
    % plot(y) 
% gdy wartoœci maj¹ kropki trzeba je zamieniæ, by móc zapisaæ w pliku
    % Estr=sprintf('%0.5f',E);
    % Estr=strrep(Estr,'.','-' );
%  nazwa(name) wykresu wa¿ne by pisaæ w "" a nie w '' bo wtdy nie dzia³a, robienie stringa w
%  zale¿noœci od parametrów wykresu
    % nazwa = "DMC_d"+D+"_N"+N+"_Nu"+Nu+"_lam"+lambda+"_e"+Estr; 
% podpis jaki bêdzie pod wykresem na overleafie
    % caption="Wykres $D=\\num{"+D+"},\\  N = \\num{"+N+"},\\  N_{\\mathrm{u}} = \\num{"+Nu+"}, \\lambda =\\num{"+lambda+"},\\  E=\\num{"+round(E,3)+"}$";
% labelka w overleafie, warto daæ ja od parametru który siê zmienia, wtedy ³atwo siê mo¿na do niej dostaæ 
    % label="fig:lambda"+lambda;
%  nazwa dokumentu w którym bêdzie skrypt do rysowania obrazków
    % fileName="lambda";
% œcie¿ka w której bêd¹ siê na overleafie rysunki znajdowaæ
    % overLeafFilePath="img/5_img/DMC/lambda/";
% zapis do pliku     
    % saveFigure(gcf, 'img/5img/DMC/',name,fileName,overLeafFilePath,caption,label);

% WA¯NA UWAGA, by nadpisaæ obecne rzeczy w pliku nale¿y odpaliæ skrypt po
% 10 sekundach od ostatniego wywo³ania, jest timer który kasuje wtedy zawartosæ
% pliku
function  saveFiguresInColumn(figures,path,names,fileName,overLeafFilePath,caption,label)
    fileContent = {  
    '\\begin{figure}[h!]'
    '   \\centering'
    };

    for i = 1:length(figures)
        saveas(figures(i),path+names(i),'epsc')

        innerFileContent={
            '   \\begin{subfigure}[b]{0.4\\textwidth}'
            '      \\includegraphics[width=1\\linewidth]{'+overLeafFilePath+names(i)+'.eps}'
            '      \\caption{}'
            '      \\label{fig:'+ label+i+'}'
            '   \\end{subfigure}'
            '       '
        };
    
%         fileContent
%         innerFileContent
        
        fileContent=[fileContent; innerFileContent];
    end

    endOFileContent={    
    '   \\caption{'+caption+'}'
    '   \\label{'+label+'}'
    '\\end{figure}'
    '           '
    };
    fileContent=[fileContent; endOFileContent];
    
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

