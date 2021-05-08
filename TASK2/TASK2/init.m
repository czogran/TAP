clc;
clear;
close all;
addpath('..\common\');

%DATA
%F-> dV/dT [cm^3/s]
%T-> temperature

%INPUTS
% Fh Fcin Fc Fd
%Hot water
Th=62;
Fh=14;
%Cold water
Tc=23;
Fcin=37;
Fc=37;
%Discruption
Td=33;
Fd=12;

Th0=62;
Fh0=14;
%Cold water
Tc0=23;
Fc0=37;
%Discruption
Td0=33;
Fd0=12;



%OUTPUT
Tout=0;
F=0;

h=81;
h0=81;

V0 = volume(h0);


T=33.57143;
T0 = T;
%ADJUSTABLE SIZES
% h;
% Tout;

%CONTROL VALUES
% Fh;
% Fcin-> Fc delayed;

%VARIABLES
% sample time
Tp=1;
t=0:Tp:10000;

%DELAY
delayT=120;
delayC=180;
ratio = 100;