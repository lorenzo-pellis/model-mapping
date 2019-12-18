function F = create_population_composition(H)

% This function generates the size distribution of a randomly selected 
% household.
% 
% Input: an (n+1)x(n+1) matrix with element (i,j)
% giving the probability that a randomly selected household has (i-1)
% adults and (j-1) children, where n is the largest household size.
% 
% Output: a 1x2 vector [F_a,F_c] for the fraction of adults and children 
% in the population.
% 
% Reference: F_a and F_c in Supplementary Methods, Section 1.5 of
% Pellis, L. et al (2019), Nature Communications
% 
% Author: Lorenzo Pellis
% Last update: 11/05/2019 

[ max_inA,max_inC ] = size(H);
F = [0,0]; % composition of the population (in fractions of total population)
F_temp = [0,0];
for inA = 1:max_inA
    nA = inA - 1;
    for inC = 1:max_inC
        nC = inC - 1;
        n = [ nA, nC ];
        F_temp = F_temp + n .* H(inA,inC);
    end
end
F = F_temp / sum(F_temp);
