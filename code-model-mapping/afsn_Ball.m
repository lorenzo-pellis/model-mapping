function afs = afsn_Ball(n,a,lambda,distribution,iota,alpha)

% This function calculates the average final size (afs), expressed in 
% numbers (hence the n of afsn) - INCLUDING the initial infectives - of a 
% within-household epidemic (or just any small population), for a single-
% type model.
% 
% Input:
%   - n: the population size
%   - a: number of initial infectives
%   - lambda: the 1-to-1 transmission rate
%   - distribution: distribution for the duration of infectious period:
%       - 0: constant, with value iota (i.e. Dirac-delta-distributed r.v.)
%       - 1: exponential distribution, with mean iota
%       - 2: Gamma distribution, with mean iota and shape parameter alpha
%   - iota: mean duration of infectious period
%   - alpha: shape parameter for the duration of infectious period (only 
%   for distribution 2)
%   - eta: exponent of frequency dependent transmission, set at 1 by
%   default in the model mapping paper, but can be changed in general
%   - den_str: a string that specifies the form of frequency dependence:
%       - 'n': transmission has the form beta/n^eta
%       - 'n-1': transmission has the form beta/(n-1)^eta
% 
% Note 1: The final size includes the initial infectives a, and the initial
% number of susceptibles is s = n - a
% 
% Note 2: The construction of the moment generating function is moved to a 
% different function, called mgf
% 
% Note 3: For efficiency, the binomial coefficients are pre-computed 
% outside this function, with nchoosek(i,j) stored in bincoeffmat(i+1,j+1)
% (to allow for 0 values for i and j).
% 
% Reference: book by Andersson & Britton (2000), Theorem 2.2
% 
% Author: Lorenzo Pellis
% Last update: 26/05/2019 

global bincoeffmat;

if n == a
    afs = a;
else
    s = n - a; % Number of initial susceptibles

    M = zeros(s+1);
    b = zeros(s+1,1);
    for h = 1:(s+1)
        for z = 1:h
            M(h,z) = bincoeffmat(s-z+2,h-z+1) / (mgf(lambda*(s-h+1),distribution,iota,alpha))^(a+z-1);
            b(h) = bincoeffmat(s+1,h);
        end
    end

    P = M \ b;
    afs = a + (0:s) * P;

    if s > 20
        warning('Warning! Population size over 20! Possible inaccuracies')
    end
end
