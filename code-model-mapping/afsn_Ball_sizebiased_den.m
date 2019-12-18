function afsn = afsn_Ball_sizebiased_den(PI,Rh,distribution,iota,alpha,eta,den_str)

% This function calculates the average final size (afs), expressed in 
% numbers (hence the n of afsn) - INCLUDING the initial infective - of a 
% within-household epidemic for a single-type model with a full household 
% size distribution. The core component is the function that computes it 
% for a single household of specified size. Then it's just a matter of 
% averaging over the correct household size distribution.
% 
% Input:
%   - PI: the size-biased household size distribution.
%   - Rh: the mean number of infectious contact an infectives
%   makes with all other household members, throughout his/her entire
%   infectious period (the within-household R0)
%   - distribution: this can be
%       - 0: for the constant iota (i.e. Dirac-delta-distributed r.v.)
%       - 1: for an exponential distribution, with mean iota
%       - 2: for a Gamma distribution, with mean iota and shape parameter alpha
%   - iota: mean of the r.v. (only for distributions 1 and 2)
%   - alpha: shape of the r.v. (only for distribution 2)
%   - eta: exponent of frequency dependent transmission, set at 1 by
%   default in the model mapping paper, but can be changed in general
%   - den_str: a string that specifies the form of frequency dependence:
%       - 'n': transmission has the form beta/n^eta
%       - 'n-1': transmission has the form beta/(n-1)^eta
% 
% Author: Lorenzo Pellis
% Last update: 26/05/2019 

max = length(PI);
temp = 0;
for n = 1:max
    if strcmp(den_str,'n-1')
        lambda = Rh / (iota*(n-1)^eta);
    else
        lambda = Rh / (iota*n^eta);
    end
    temp = temp + PI(n) * afsn_Ball(n,1,lambda,distribution,iota,alpha);
end
afsn = temp;

