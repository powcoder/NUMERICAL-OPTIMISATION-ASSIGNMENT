function [d2f] = computeD2f(ts, yst, a, b, c )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

residual = model(ts, a, b, c) - yst;

jaco = computeJacobian(ts, a, b, c);

d2f = jaco' * jaco;

m = length(ts);

for i=1:m
    
    t = ts(i);
    
    z = [0  						0              -t * exp(-c*t)
            0  						0              -t^3 * exp(-c*t)
        t*(-exp(-c*t))           -t^3 * exp(-c*t)   t^2 * exp(-c*t)* (a + b*t^2)];
    
    d2f = d2f + residual(i) * z;
end

end

