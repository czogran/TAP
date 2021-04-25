%% GPC
clc
close all
clear all

sympref('FloatingPointOutput',true)
syms uMin [1,15]
syms yMin [1,15]
syms y0 
syms y [1,15]
syms yZad

a=0.2676;
% INDEXES FROM 1-> b0=b1
b1=0;
b2=0;
b3=0.1989;
b4=0.2552;
b5=0;
b6=0;

startPoint=10;

N=4;
Nu=2;
lambda=0.1;
psi=1;

yVector=cell(2*startPoint,1);
yVector{startPoint-1}=yMin(1);
yVector{startPoint}=y0;

d=yVector{startPoint}-a*yVector{startPoint-1}-...
    b1*uMin(1)-...
    b2*uMin(2)-...
    b3*uMin(3)-...
    b4*uMin(4)-...
    b5*uMin(5)-...
    b6*uMin(6);

index=1;
for k=startPoint:startPoint+N

    yVector{k+1}=a*yVector{k}+...
        b1*uMin(max(1,1-index))+...
        b2*uMin(max(1,2-index))+...
        b3*uMin(max(1,3-index))+...
        b4*uMin(max(1,4-index))+...
        b5*uMin(max(1,5-index))+...
        b6*uMin(max(1,6-index))+d;
    
    index = index + 1;
end

stepResponses={N+1};
for  k=1:N+1
    stepResponses{k} = yVector{startPoint+k-1};
end

stepResponsesN={N};
for  k=1:N
    stepResponsesN{k} = yVector{startPoint+k};
end

% SUBSTITUTIONS
% ACHRUNG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% VERY IMPORTANT TO SET 0 FOR DELAYD VALUES
uMin1=1;    yMin1=0; 
uMin2=0;    yMin2=0; 
uMin3=0;    yMin3=0; 
uMin4=0;    yMin4=0; 
uMin5=0;    yMin5=0; 
uMin6=0;    yMin6=0; 
uMin7=0;    yMin7=0; 
uMin8=0;    yMin8=0; 
uMin9=0;    yMin9=0; 
uMin10=0;    yMin10=0; 
uMin11=0;    yMin11=0; 
uMin12=0;    yMin12=0; 
uMin13=0;    yMin13=0; 
uMin14=0;    yMin14=0; 
uMin15=0;    yMin15=0; 

y0=0;
y=0;

% EVALUATE
steps= [N];
for k=1:N+1
    steps(k) = subs(stepResponses{k})
end

M=MMatrix(N,Nu,1,1,steps)

psiMatrix=diag(ones(N,1))*psi;
lambdaMatrix=diag(ones(Nu,1))*lambda;

K=(M'*psiMatrix*M+lambdaMatrix)^(-1)*M'*psiMatrix

K1=K(1,:);


law = 0;
for i =1:N
    law= law+K1(i)*(yZad-stepResponsesN{i});
end
law