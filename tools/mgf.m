function mgf = mgf(t,distribution,iota,alpha)

% This function computes the moment-generating function, defined for r.v. X
% as mgf(t)=E[exp(-tX)], as probabilists tend to define it (i.e. like a
% Laplace transform) for 3 main distributions used in the standard
% stochastic SIR model.
% 
% Input:
%   - t: the argument of the mgf
%   - distribution: this can be
%       - 0: for the constant iota (i.e. Dirac-delta-distributed r.v.)
%       - 1: for an exponential distribution, with mean iota
%       - 2: for a Gamma distribution, with mean iota and shape parameter alpha
%   - iota: mean of the r.v. (only for distributions 1 and 2)
%   - alpha: shape of the r.v. (only for distribution 2)
% 
% Author: Lorenzo Pellis
% Last update: 18-05-2019 

if distribution==0
    mgf=exp(-t*iota);
elseif distribution==1
        mgf=1/(1+iota*t);
else
        mgf=(alpha/(alpha+iota*t))^alpha;
end

