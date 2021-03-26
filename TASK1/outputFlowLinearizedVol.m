function F=outputFlowLinearizedVol(V,V0)
    %a-> constant
    a=7;
    C=0.7;
    F=a*(sqrt(sqrt(V0/C))+(0.25/sqrt(sqrt(C*(V0^3))))*(V-V0));
end

