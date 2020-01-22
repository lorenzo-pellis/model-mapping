function esc_prob = esc_prob_multitype(h,Lambda,z)

% This function computes the probability that a randomly selected
% susceptible escapes global infection from all infectious individuals,
% expressed as fractions of each class by the components of vector z. This
% is relevant to compute the asymptotic final size for a multitype model,
% both with and without households.
%
% Input:
%   - h: vector with elements given the fractions of individuals of each
%   type in the population (if total is not 1, it is normalised)
%   - Lambda: matrix with the total infectivities towards all individuals
%   in the population (i.e. 1-to-all infectivities, given by 1-to-all
%   infection rates times average duration of infection). Therefore,
%   element l_{ij} is such that l_{ij}/N is the 1-to-1 infectivity, i.e.
%   the total infectivity an individual of type j exerts on a single
%   susceptible of type i, throughout the entire infectious period. Note
%   that if for a multitype households model, Lambda should ouly contain
%   global infectivities.
%   - z: column vector giving the fraction of individuals of each type that
%   are infectious. These are not fractions of the total population, but of
%   each group separately (i.e. z_i = Z_i/N_i, where N_i is the size of 
%   group i)
% 
% Note 1: Lambda is not the NGM, but rather NGM = ( h_i l_{ij} ). I could
% have passed the NGM, but in this was the computation of Lambda is
% offloaded outside the function (which is called multiple times when
% solving the final size implicit equations). I actually tried to do the
% computations all in one line (esc_prob = exp(-((1./h).*(NGM*(h.*z))))),
% but it is actually slower than the double for loop below (at least in
% dimension 2). There is a different version of this function, called 
% "esc_prob_multitype_NGM", that takes the NGM as input and performs the 
% computation above.
% 
% Note 2: In the references below, Lambda is typically transposed to what
% is used here. Here the first index is the infectee and the second the
% infector.
% 
% References:
% For multitype model: book by Andersson & Britton (2000), Section 6.2
% For multitype-households models: Ball, Britton & Sirl (2011), Journal of 
% Mathematical Biology
% 
% Author: Lorenzo Pellis
% Last update: 26-05-2019 

ntypes = length(h); % Number of types
[check1,check2] = size(Lambda);
check3 = length(z);
if check1 ~= check2
    error('Error: Next generation matrix non squared!');
end
if check1 ~= ntypes
    error('Error: number of types does not agree with the size of the next generation matrix!');
end
if check3 ~= ntypes
    error('Error: number of types does not agree between vectors!');
end
if ~(sum(h)==1)
    warning('Warning: population group fractions don''t sum to 1, so I''m renormalising...')
    h = h / sum(h); % Renormalise into fractions of total population
end

esc_prob = NaN(ntypes,1);
for n = 1:ntypes
    if h(n)==0
        esc_prob(n,1) = 1;
    else
        temp = 0;
        for m = 1:ntypes
            temp = temp + Lambda(n,m)*h(m)*z(m);
        end
        esc_prob(n,1) = exp(-temp);
    end
end
