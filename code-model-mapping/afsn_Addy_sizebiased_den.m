function afs = afsn_Addy_sizebiased_den(PI,pesc,Rh,distribution,iota,alpha,eta,den_str)

% This function calculates the average final size (afs), expressed in 
% numbers (hence the n of afsn) of an epidemic in the household of a 
% randomly selected individual, for a single-type model with a full 
% household size distribution, when individuals are allowed to escape 
% infection independently from outside and there are no initial infectives.
% The core component is the function that computes it for a single 
% household of specified size. Then it's just a matter of averaging over 
% the correct (i.e. size-biased) household size distribution.
% 
% Input:
%   - PI: the size-biased household size distribution (single-type)
%   - pesc: probability that each single individual escapes infection from
%   outside
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
% Output: average final size (numbers)
% 
% Note 1: For efficiency, the binomial coefficients are pre-computed 
% outside this function, with nchoosek(i,j) stored in bincoeffmat(i+1,j+1)
% (to allow for 0 values for i and j).
% 
% Reference: Supplementary Methods, Section 1.4 of
% Pellis, L. et al (2020), Nature Communications
%
% Methodological references: 
% Addy, Longini and Haber (1991), Biometrics
% Ball and Lyne (2001), Advances in Applied Probability
% 
% Author: Lorenzo Pellis
% Last update: 14-06-2019 

nmax = length(PI);
temp = 0;
for n = 1:nmax
    if strcmp(den_str,'n-1')
        if n == 1
            lambda = 0; % It doesn't matter, as within-household transmission is a household of size 1 is irrelevant
        else
            lambda = Rh/(n-1)^eta;
        end
    else
        lambda = Rh/n^eta;
    end
    temp = temp + PI(n) * afsn_Addy(n,pesc,lambda,distribution,iota,alpha) / n;
end
afs = temp;

end
