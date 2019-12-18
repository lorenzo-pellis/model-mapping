function afsn = afsn_Ball_incl(n,a,lambda,distribution,iota,alpha)

% This function calculates the average final size (afs), expressed in 
% numbers (hence the n of afsn) - INCLUDING the initial infectives - of a 
% within-household epidemic (or just any small population), for a single-
% type model. There is an identical function called "afsn_Ball" that uses a
% slightly different system, where the auxiliatory variable runs from 0 to
% s (initial number of susceptibles) rather than a (initial nuber of
% infectives) to n (total group size).
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

if n == a
    afsn = a;
else
    global bincoeffmat; %#ok<TLEV>

    A = zeros(n-a+1);
    b = zeros(n-a+1,1);

    ia = a+1;
    for il = ia:(n+1)
        l = il-1;
        for iu = ia:il
            u = iu-1;
            A(il-a,iu-a) = bincoeffmat(n-u+1,l-u+1) / (mgf(lambda*(n-l),distribution,iota,alpha))^u;
        end
        b(il-a) = bincoeffmat(n-a+1,l-a+1);
    end
    P = A\b;

    afsn = (a:n) * P;

end