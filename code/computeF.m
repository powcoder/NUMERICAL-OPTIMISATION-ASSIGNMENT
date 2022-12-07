function [ f ] = computeF(ts, yst, a, b, c)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

residual = model(ts, a, b, c) - yst;

f =  0.5 * sum(residual.^2);

end

