function y = draw_random_deviates_from_empirical_cum_distribution(cumdistr,dim,support)

% cumdistr is a column vector of containing the cumulative distribution of
% the empirical discrete distribution to draw random deviates from. It is
% assumed that the first value is the probability of drawing a random
% deviate from support(1) and that cumdistr(end) = 1. It is assumed the
% values in cumdistr are in increasing order.

% If cumdistr(end)<1, a warning is given and NaNs are generated with the
% excess probability.

% If support is not provided, {1,2,...,length(cumdistr)} is the assumed
% support

if ~isvector(cumdistr) || ( cumdistr(end) > 1 ) || any( cumdistr ~= sort(cumdistr) )
    disp('Abort! Invalid cumulative distribution:')
    disp(cumdistr);
    beep;
elseif cumdistr(end) < 1
    disp(['Defective cumulative distribution. NaNs are generated with probability ',num2str(1-cumdistr(end))]);
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

