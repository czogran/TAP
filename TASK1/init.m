clc;
clear;
close all;
addpath('..\common\');

% SETTINGS
colorLabels=['m','c','r','g','b','w','k','y'];

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
%Discruption
Td=33;
Fd=12;

%OUTPUT
Tout=0;
F=0;

h=81;
h0=81;
T=33.57;
T0=33.57;

%ADJUSTABLE SIZES
% h;
% Tout;

%CONTROL VALUES
% Fh;
% Fcin-> Fc delayed;

%VARIABLES
% sample time
Tp=1;
t=0:Tp:2000;

%DELAY
delayT=120;
delayC=180;