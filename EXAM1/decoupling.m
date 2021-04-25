%% DECOUPLING
clc
clear all
close all

syms s

G11= tf([1],[24 10 1]);
G11Dealay=exp(-5*s);

G12 =tf([0.4],[9 1]);
G12Delay=exp(-7*s);

G21=tf([0.5],[5 1]);
G21Dealay=exp(-6*s);

G22=tf([-0.4 2],[4 1]);
G22Dealay=exp(-7*s);


D12=-G12/G11;
D12Delay= G12Delay/G11Dealay;
D12DelayClone =D12Delay;

D21=-G21/G22;
D21Delay= G21Dealay/G22Dealay;
D11DelayClone =D21Delay;

disp("--------------------------------------------------------------------")
disp("D12=-G12/G11");
D12
s = 1;
if length(zero(D12))>length(pole(D12))
    D12MessageRate=("DO MIANOWNIKA D12 Dodajemy inercję, o odpowiednio małej stałej czasowej (o rząd wielkości mniejszej od najmniejszej stałej czasowej)")
end

if find(pole(D12)<2)
    D21MessageRate=("Z MIANOWNIKA D21 wywalamy ujemne")
end

if subs(D12Delay)> 1
    D12MessageDelay=("Z D12 trzeba wywalić stałą czasową, poniewaz przewiduje w przód")
end

disp("--------------------------------------------------------------------")
disp("--------------------------------------------------------------------")
disp("D21=-G21/G22")
D21

if length(zero(D21))>length(pole(D21))
    D21MessageRate=("DO MIANOWNIKA D21 Dodajemy inercję, o odpowiednio małej stałej czasowej (o rząd wielkości mniejszej od najmniejszej stałej czasowej)")
end

if find(pole(D21)<2)
    D21MessageRate=("Z MIANOWNIKA D21 wywalamy ujemne")
end

if subs(D21Delay)> 1
    D21MessageDelay=("Z D21 trzeba wywalić stałą czasową, poniewaz przewiduje w przód")
end