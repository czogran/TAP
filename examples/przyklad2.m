%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                              %
%%    Skrypt ilustruj¹cy linearyzacjê modelu zbiornika                          %
%%  (przyk³ad II z wprowadzenia).                                               %
%%    h - wysokoœæ s³upa cieczy w zbiorniku                                     %
%%    hL - wyjœcie zlinearyzowanego modelu zbiornika                            %
%%    F - dop³yw cieczy do zbiornika (mo¿na zmieniaæ w linii 36.)               %
%%                                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Autor: Piotr Marusak; 2002-2021
%   Wersja przeznaczona dla s³uchaczy wyk³adu Technika Automatyzacji Procesów.

clear all;

% parametry zbiornika

P=530;
alfa=27;

% warunki pocz¹tkowe

F=111;
h=(F/alfa)^2;
h0=h;
hL=h;
V=P*h;
VL=P*hL;

% okres próbkowania

Tp=1;

% skokowe zmiany sterowania

F=F+10; %%%%% warto pozmieniaæ (np. o -30, -20, -10, 10, 20, 30)

for i=1:2000
   
   % rozwi¹zywanie równañ ró¿niczkowych zmodyfikowan¹ metod¹ Eulera
   
   pomV=V+0.5*Tp*(F-(alfa*sqrt(h)));
   pomh=pomV/P;
   
   dV=Tp*(F-(alfa*sqrt(pomh)));
   V=V+dV;
   
   pomVL=VL+0.5*Tp*(F-(alfa*(sqrt(h0)+(1/(2*sqrt(h0)))*(hL-h0))));
   pomhL=pomVL/P;
   
   dVL=Tp*(F-(alfa*(sqrt(h0)+(1/(2*sqrt(h0)))*(pomhL-h0))));      
   VL=VL+dVL;
   
   % przygotowanie zmiennych do rysowania
   
   wyV(i)=V;
   wyVL(i)=VL;
   
   h=V/P;
   hL=VL/P;
   wyh(i)=h;
   wyhL(i)=hL;
   
   czas(i)=(i-1)*Tp;
end

% rysowanie przebiegów zmiennych

plot(czas,wyh,'b');hold on;plot(czas,wyhL,'m');
xlabel('czas');
ylabel('h');
