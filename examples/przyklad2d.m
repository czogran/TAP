%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                              %
%%    Skrypt ilustruj¹cy dyskretyzacjê transmitancji opisuj¹cej zbiornik        %
%%  (ci¹g dalszy przyk³adu II z wprowadzenia; najpierw nale¿y uruchomiæ         %
%%  przyklad2.m).                                                               %
%%    Tp - stopniowo zmieniany okres próbkowania                                %
%%         (kolejne zmiany po naciœniêciu dowolnego klawisza)                   %
%%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Autor: Piotr Marusak; 2013-2021
%   Wersja przeznaczona dla s³uchaczy wyk³adu Technika Automatyzacji Procesów.

Tp=200; % pocz¹tkowa wartoœæ okresu próbkowania (wyraŸnie zbyt du¿a)
metoda='zoh';
sys=tf(1,[P alfa/(2*sqrt(h0))]) % transmitancja ci¹g³a
sysd=c2d(sys,Tp,metoda) % transmitancja dyskretna
step(sys); hold on; step(sysd)

pause;
Tp=Tp/2
sysd1=c2d(sys,Tp,metoda); step(sysd1);

pause;
Tp=Tp/2
sysd2=c2d(sys,Tp,metoda); step(sysd2);

pause;
Tp=Tp/2.5
sysd3=c2d(sys,Tp,metoda); step(sysd3);

pause;
Tp=Tp/2
sysd4=c2d(sys,Tp,metoda); step(sysd4);

pause;
Tp=Tp/10
sysd5=c2d(sys,Tp,metoda); step(sysd5);
