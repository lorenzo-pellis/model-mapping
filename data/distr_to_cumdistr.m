function y = distr_to_cumdistr(x)

% This is a handmade function that takes a 1D empirical probability
% mass function (like a histogram count) and creates a cumulative mass 
% function for the same distribution.
%   
% Author: Lorenzo Pellis
% Last update: 21-01-2020 

if ~isvector(x)
    disp('Abort! Invalid distribution:');
    disp(x);
    beep;
    return;
end

n = length(x);
y = zeros(size(x));
y(1) = x(1);
for i = 2:n
    y(i) = y(i-1) + x(i);
end

