function F=outputFlowLinearized(h,h0)
    %a-> constant
    a=7;
    F=a*(sqrt(h0)+(0.5/sqrt(h0))*(h-h0));
end


