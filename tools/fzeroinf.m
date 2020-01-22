function z = fzeroinf(fun,lb,options,factor,maxiter,initguess)

% This is my own variant of fzero, to be used when I know there is a
% lower bound on the solution, but where the function is not defined:
% lb = lower bound, which is an inf here, rather than a min
% 
% Note that initguess is not meant to be a likely upper bound for the
% solution. Rather, I start searching from initguess and gradually expand
% the interval squashing towards lb on the left and running to infinity on
% the right.
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

if nargin < 5 || isempty(maxiter)
    maxiter = 100;
end
if nargin < 4 || isempty(factor)
    factor = 10;
end
if nargin < 6;
    il = lb + 1;
else
    il = initguess;
end

ir = il;
signl = sign(fun(il));
signr = sign(fun(ir));
% disp([0,fun(il),fun(ir)])
if signr ~= signl
    if ( nargin < 3 || isempty(options) )
        z = fzero(fun,[il,ir]);
    else
        z = fzero(fun,[il,ir],options);
    end
    return;
else
    flag = false;
    for n = 1:maxiter
        l = lb + (il-lb)/factor^n;
        r = ir + factor^n;
        signl = sign(fun(l));
        signr = sign(fun(r));
%         disp([n,fun(il),fun(ir)])
        if signr ~= signl
            flag = true;
            break;
        end
    end
    if flag
        if ( nargin < 3 || isempty(options) )
            z = fzero(fun,[l,r]);
        else
            z = fzero(fun,[l,r],options);
        end
    else
        lx = zeros(1,maxiter);
        ly = zeros(1,maxiter);
        rx = zeros(1,maxiter);
        ry = zeros(1,maxiter);
        for n = 1:maxiter
            l = lb + (il-lb)/factor^n;
            lx(n) = l;
            ly(n) = fun(l);
            r = ir + factor^n;
            rx(n) = r;
            ry(n) = fun(r);
        end
        z = NaN;
        disp('fzeroinf failed! Check your function...');
        disp('Intermediate calculations:');
        disp([il,lx]);
        disp([ir,rx]);
        disp([fun(il),ly]);
        disp([fun(ir),ry]);
    end
end
    