function [ jaco ] = computeJacobian(ts, a, b, c)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

m = length(ts);

jaco = zeros(m, 3);

for i=1:m
    t = ts(i);
    jaco(i, 1) = exp(-c*t);
    jaco(i, 2) = t^2 * exp(-c*t);
    jaco(i, 3) = - t * exp(-c*t) * (a + b* t^2);
    
end


end

