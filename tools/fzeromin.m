function z = fzeromin(fun,lb,options,factor,maxiter,initguess)

% This is my own variant of fzero, to be used when I know there is a
% minimum value of the solution, on which the function is defined:
% lb = lower bound, which is also a min here
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
% Last update: 12/05/2019 

flag = false;
if nargin < 5 || isempty(maxiter)
    maxiter = 100;
end
if nargin < 4 || isempty(factor)
    factor = 10;
end
signlb = sign(fun(lb)); % The function is defined at the lower bound
if nargin < 6
    l = lb;
    signl = signlb;
else
    signi = sign(fun(initguess));
    if signi ~= signlb % Found the starting interval already
        if ( nargin < 3 || isempty(options) )
            z = fzero(fun,[lb,initguess]);
        else
            z = fzero(fun,[lb,initguess],options);
        end
        return
    else
        l = initguess;
        signl = signi;
    end
end

% Procedure to find the upper bound of the interval...
r = l;
for n = 1:maxiter
    r = r + factor^n;
    signr = sign(fun(r));
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
    rvec = zeros(maxiter,2);
    for n = 1:maxiter
        rvec(n,1) = l + factor^n;
        rvec(n,2) = fun(rvec(n,1));
    end
    z = NaN;
    disp('fzeroup failed! Check your function...');
    disp('Intermediate calculations:');
    disp([lb,fun(lb)]);
    disp([l,fun(l)]);
    disp(rvec);
end
    