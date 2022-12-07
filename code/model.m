function [y] = model(t, x1, x2, x3)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

y = (x1 + x2* (t.^2)) .* exp(-x3 * t);

end





