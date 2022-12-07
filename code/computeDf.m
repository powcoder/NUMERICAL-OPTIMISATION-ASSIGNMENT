function [ df ] = computeDf(ts,yst,  a, b, c)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

residual = model(ts, a, b, c) - yst;

jaco = computeJacobian(ts, a, b, c);

df = jaco' * residual';

end

