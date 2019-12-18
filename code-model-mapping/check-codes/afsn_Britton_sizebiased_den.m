function afs = afsn_Britton_sizebiased_den(PI,pesc,Rh,distribution,iota,alpha,eta,den_str)

% This function calculates the average final size (afs), expressed in 
% numbers (hence the n of afsn) of an epidemic in the household of a 
% randomly selected individual, for a single-type model with a full 
% household size distribution, when individuals are allowed to escape 
% infection independently from outside and there are no initial infectives.
% 
% This function is used to calculate the asymptotic final size for a
% households model with a household size distribution, and relies of the
% method described in the book by Andersson and Britton (2000), which is
% based on computing the within-household final size distribution with "a"
% initial infectives, which are those that escape infection from outside.
% This method is less efficient than finding directly the final size
% distribution with individuals escaping infection from outside and no
% initial infectives (using Addy, Longini and Haber (1991) or Ball and Lyne
% (2001)). Therefore I only use this function to cross-check the results
% more efficiently computed with the function "afsH_varHsize_Addy_den".
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
% Methodological reference: book by Andersson & Britton (2000), Section 6.3
% 
% Author: Lorenzo Pellis
% Last update: 14/06/2019 

global bincoeffmat;

nmax = length(PI);
temp = 0;
for n = 1:nmax
    if strcmp(den_str,'n-1')
        lambda = Rh/(n-1)^eta;
    else
        lambda = Rh/n^eta;
    end
    temp2 = 0;
    for a = 1:n
        temp2 = temp2 + bincoeffmat(n+1,a+1) * pesc^(n-a) * (1-pesc)^a * afsn_Ball(n,a,lambda,distribution,iota,alpha) / n;
    end    
    temp = temp + PI(n) * temp2;
end
afs = temp;

end
