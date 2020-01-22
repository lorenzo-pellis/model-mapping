function y = draw_random_deviates_from_empirical_cum_distribution(cumdistr,dim,support)

% This function draws random numbers from an empirical cumulative
% mass function (discrete).
% 
% Inputs:
%   - cumdistr: a column vector of containing the cumulative distribution 
%       of the empirical discrete distribution to draw random deviates 
%       from. It is assumed that the first value is the probability of 
%       drawing a random deviate from support(1) and that cumdistr(end)=1. 
%       It is assumed the values in cumdistr are in increasing order.
%   - dim: number of dimensions of the distribution. If not provided, it is
%       assumed to be 1.
%   - support: the values for which the elements of cumdistr give the
%       cumulative mass function. If not provided, the support is assumed
%       to be {1,2,...,length(cumdistr)}.
% 
% Notes: If cumdistr(end)~=1, a warning is given and, if defective, NaNs 
% are generated with the missing probability.
% 
% Update: 21-01-2020

if ~isvector(cumdistr)
    disp(cumdistr);
    error('Abort! Invalid cumulative distribution:')
end
% In case the cumulative distribution might not be sorted, uncomment the
% line below:
% cumdistr = sort(cumdistr)
if cumdistr(end) > 1
    warning(['Cumulative distribution is excessive by a mass of ', num2str(cumdistr(end)-1), ' and the excess mass is lost!'])
elseif cumdistr(end) < 1
    warning(['Defective cumulative distribution. NaNs are generated with probability ',num2str(1-cumdistr(end))]);
end

n = length(cumdistr);
if nargin < 2
    dim = 1;
    support = 1:n;
elseif nargin < 3
    support = 1:n;
end

D = prod(dim);
r = rand(D);
index = zeros(D);
for i = 1:D
    index(i) = find(cumdistr >= r(i),1);
end
if ( length(dim) > 1 ) && ( dim(2) > 1 )
    y = reshape(support(index),dim);
else
    y = support(index);
end

