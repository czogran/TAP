function h=heightFromVolume(V)
    %C-> constant
    C=0.7;    
    if(V==0)
        h=0;
    else
        h=sqrt(V/C);
    end
end