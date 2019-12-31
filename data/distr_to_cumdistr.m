function y = distr_to_cumdistr(x)

% Assumes a 1D discrete distribution (like a histogram count)

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

