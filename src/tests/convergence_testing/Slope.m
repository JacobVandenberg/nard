function [m] = Slope(x, y)
%SLOPE Summary of this function goes here
%   Detailed explanation goes here
    F = [x.'.^0 x.'];
    c = F\ (y.');
    m = c(2);
end

