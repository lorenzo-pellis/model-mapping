function R0 = R0H(chain,Rg)

% This function computes R0 for a model with households.
%
% Input:
%   - chain: this is a vector of length N, giving the average number of 
%   cases in each generation. The average is over the stochastic dynamics 
%   and the household size distribution, so N is the largest possible 
%   number of generations of all household size allowed by the population 
%   structure.
%   - Rg is the average number cases a single case makes through global 
%   contacts only, during the entire infectious period, in a fully 
%   susceptible population, i.e. the basic reproduction number for global 
%   contacts only.
%
% Note 1: R0 is computed as the dominant eigenvalue of a matrix, rather
% than by solving the characteristic equation directly.
% 
% Reference: Supplementary Methods, Section 1.4 of
% Pellis, L et al (2020), Nature Communications
% 
% Methodological reference: Pellis, Ball Trapman (2012), Mathematical
% Biosciences, Theorem 1.
% 
% Author: Lorenzo Pellis
% Last update: 29-12-2019 

nmax = length(chain);
M = zeros(nmax);
M(1,:) = Rg*chain;
for i = 2:nmax
    M(i,i-1) = 1;
end

R0 = max(eig(M));

