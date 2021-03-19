%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                              %
%%    Skrypt ilustruj¹cy linearyzacjê przyk³adowej funkcji                      %
%%  (przyk³ad Ia z wprowadzenia).                                               %
%%    x0 - punkt wzglêdem którego jest wykonywana linearyzacja (linia 14.)      %
%%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Autor: Piotr Marusak; 2013-2021
%   Wersja przeznaczona dla s³uchaczy wyk³adu Technika Automatyzacji Procesów.

clear all;

x0=-1 %%%%% warto pozmieniaæ
x=-2:0.01:2;
y=2.*(x.^2)+4;
yL=2.*x0.^2+4+4.*x0.*(x-x0);

plot(x,y,'b'); hold on; xlabel('x'); ylabel('y'); % rysunek funkcji nieliniowej
plot(x0,2.*(x0.^2)+4,'og') % punkt, w którym funkcja by³a linearyzowana
plot(x,yL,'g'); % rysunek funkcji liniowej
