function y = normalise_distribution(x)

if ~isvector(x)
    disp('Abort! Invalid distribution:');
    disp(x);
    beep;
    return;
end

y = x/sum(x);