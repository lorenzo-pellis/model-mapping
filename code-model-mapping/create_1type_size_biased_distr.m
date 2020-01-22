function PI_single = create_1type_size_biased_distr(H_single)

% This function generates the size-biased household size distribution, 
% i.e. the size distribution of the household of a randomly selected 
% individual (irrespective of age class).
% 
% Input: 1xn vector of the size distribution of a randomly selected 
% household.
% 
% Output: a 1xn vector, with i-th element giving the probability that 
% a randomly selected household has size i.
% 
% Reference: \pi_n in Supplementary Tables 2, 5 and 7 of 
% Pellis, L. et al (2020), Nature Communications
% 
% Author: Lorenzo Pellis
% Last update: 11-05-2019 


max_size = length(H_single);
PI_single = zeros(1,max_size);
PI_single_temp = PI_single;

for n = 1:max_size
    PI_single_temp(n) = n * H_single(n);
end
PI_single = PI_single_temp / sum(PI_single_temp);
