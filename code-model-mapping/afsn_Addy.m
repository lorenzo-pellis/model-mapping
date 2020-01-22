function afsn = afsn_Addy(n,pesc,lambda,distribution,iota,alpha)

% This function calculates the average final size (afs), expressed in 
% numbers (hence the n of afsn) of an epidemic in a small population (e.g. 
% a household, for a single-type model where individuals escape infection 
% from outside independently with specified probability and there are no
% infectives at the start.
% 
% Input:
%   - n: the size of the group (all individuals initially susceptible)
%   - pesc: probability that each single individual escapes infection from
%   outside
%   - lambda: the 1-to-1 transmission rate
%   - distribution: distribution for the duration of infectious period:
%       - 0: constant, with value iota (i.e. Dirac-delta-distributed r.v.)
%       - 1: exponential distribution, with mean iota
%       - 2: Gamma distribution, with mean iota and shape parameter alpha
%   - iota: mean duration of infectious period
%   - alpha: shape parameter for the duration of infectious period (only 
%   for distribution 2)
% 
% Output: average final size (numbers)
% 
% Note 1: For efficiency, the binomial coefficients are pre-computed 
% outside this function, with nchoosek(i,j) stored in bincoeffmat(i+1,j+1)
% (to allow for 0 values for i and j).
% 
% Note 2: For generality, any distribution for the duration of the 
% infectious period is allowed. The moment-generating function of the 
% random duration is relegated to an external function "mgf"
% 
% Reference: Supplementary Methods, Section 1.4 of
% Pellis, L et al (2020), Nature Communications
%
% Methodological references: 
% Addy, Longini and Haber (1991), Biometrics
%
% Author: Lorenzo Pellis
% Last update: 14-06-2019 


global bincoeffmat;

A = zeros(n+1);
b = zeros(n+1,1);

for il = 1:(n+1)
    l = il-1;
    for iu = 1:il
        u = iu-1;
        A(il,iu) = bincoeffmat(n-u+1,l-u+1) / pesc^(n-l) / (mgf(lambda*(n-l),distribution,iota,alpha))^u;
    end
    b(il) = bincoeffmat(n+1,l+1);
end
P = A\b;

afsn = (0:n) * P;

