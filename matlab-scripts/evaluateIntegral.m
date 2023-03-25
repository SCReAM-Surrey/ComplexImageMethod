function [output] = evaluateIntegral(c,B,t)
% evaluate integral of spherical wave at a certain distance from reflector
% given limits
    
%     %form with incomplete gamma function - blows up for complex impedance
    output = exp(-1i.*B.*c).*sqrt(c.*(t - 1i*B).*igamma(0.5,c.*(t-1i*B)))./...
        (sqrt(B + 1i*t));

    %alternate form with exponential integral
%     addpath('genexpint/.');
%     output = 1i*exp(-1i*B.*c).*sqrt(B+1i*t).*genexpint(0.5,c.*(t-1i*B));
end
