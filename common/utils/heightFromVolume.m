function h=heightFromVolume(V)
    %C-> constant
    C=0.7;    
    h=sqrt(V/C);
end

%     h=sqrt(V0/C)+1*sqrt(1/C)/(2*sqrt(V0))*(V-V0);
