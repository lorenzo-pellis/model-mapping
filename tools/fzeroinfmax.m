function z = fzeroinfmax(fun,lb,rb,options,factor,maxiter,initguess)

% This is my own variant of fzero, to be used when I know the solution
% belongs to an interval, but not both extremes can (or should) be taken.
% In this case:
%   - lb (= left bound) is an inf (the function might not be defined there
%   or I might not want to calculate it there - e.g. final size equation 
%   has always 0 as a solution but I want the strictly positive one)
%   - rb (= right bound) is a max, so the function is defined there
% 
% This is because sometimes the search when fzero starting with a single
% initial guess might go wrong, so it's useful to automatically search for
% a good interval where the function changes sign. If this fails, at least
% I get an idea of which points this code has looked at and acquire some
% knowledge of the function fun and why my search fails.
% 
% All variants involve:
%   - fzeromin: when there is a lower bound where the function is defined
%   - fzeromax: when there is an upper bound where the function is defined
%   - fzeroinf: with lower bound, but where the function is not defined
%   - fzerosup: with upper bound, but where the function is not defined
%   - fzerosminsup: with two bounds and function defined only on lower
%   - fzerosinfmax: with two bounds and function defined only on upper
%   - fzerosinfsup: with two bounds and function not defined on either
% When both min and max exist, just use fzero...
% 
% Author: Lorenzo Pellis
% Last update: 06-11-2019 

if nargin < 6 || isempty(maxiter)
    maxiter = 100;
end
if nargin < 5 || isempty(factor)
    factor = 10;
end
signrb = sign(fun(rb));
if nargin < 7;
    ir = rb;
    signr = signrb;
else
    signi = sign(fun(initguess));
    if signi ~= signrb % Found the starting interval already
        if ( nargin < 4 || isempty(options) )
            z = fzero(fun,[initguess,rb]);
        else
            z = fzero(fun,[initguess,rb],options);
        end
        return
    else
        ir = initguess;
        signr = signi;
    end
end

il = ir;
r = ir;
flag = false;
for n = 1:maxiter
    l = lb + (il-lb)/factor^n;
    signl = sign(fun(l));
    if signr ~= signl
        flag = true;
        break;
    end
end
if flag
    if ( nargin < 4 || isempty(options) )
        z = fzero(fun,[l,r]);
    else
        z = fzero(fun,[l,r],options);
    end
else
    lx = zeros(1,maxiter);
    ly = zeros(1,maxiter);
    for n = 1:maxiter
        l = lb + (il-lb)/factor^n;
        lx(n) = l;
        ly(n) = fun(l);
    end
    z = NaN;
    disp('This routine failed! Check your function...');
    disp('Intermediate calculations:');
    disp([il,lx]);
    disp([fun(il),ly]);
end
    