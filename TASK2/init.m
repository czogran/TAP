clc;
clear;
close all;
addpath('..\common\');

% SETTINGS
colorLabels=['m','c','r','g','b','k','y'];

%DATA
%F-> dV/dT [cm^3/s]
%T-> temperature

%INPUTS
% Fh Fcin Fc Fd
%Hot water
Th=62;
Th0=62;

Fh=14;
Fh0=14;
%Cold water
Tc=23;
Tc0=23

Fcin=37;
Fcin0=37;
Fc0=37;
%Discruption
Td=33;
Td0=33;
Fd=12;
Fd0=12;
%OUTPUT
Tout=0;
F=0;

h=81;
h0=81;
T=33.57;
T0=33.57;

V0=0.7^h0^2;

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