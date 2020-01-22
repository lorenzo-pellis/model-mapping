function PI_type = create_type_biased_distr(PI)

% This function generates the age class-biased household size distribution, 
% i.e. the size distribution of the household a randomly selected adult/child.
% 
% Input: a (n+1)x(n+1)x2 double matrix. The first layer has elements (i,j)
% giving the probability that the household of a randomly selected adult
% has (i-1) adults and (j-1) children. The second layer for the household
% of a randomly selected child. Note: row 1 in layer 1 is 0 and column 1 in
% layer 2 is 0.
% 
% Output: a 2xn matrix. Element i in row 1 gives the probability that the
% household of a randomly selected adult has size i. Row 2 for a randomly
% selected child (check: element (2,1) should be very small or 0, if
% children are unlikely to live on their own)
% 
% Reference: \pi_n^a and \pi_n^c from Supplementaty Tables 2, 5 and 7 of
% Pellis, L. et al (2020), Nature Communications
% 
% Author: Lorenzo Pellis
% Last update: 11-05-2019
% This function creates the average size of a household of a randomly
% selected adult or child: mu = [ mu_a, mu_c ]

[ max_inA,max_inC,~ ] = size(PI);
max_size = max_inA-1; % Households are truncated at a maximum size, not max numbers of adults and children
PI_type_temp = zeros(2,max_size);
for type = 1:2
    for inA = 1:max_inA
        nA = inA - 1;
        for inC = 1:max_inC
            nC = inC - 1;
            nH = nA + nC;
            if ( nH ~= 0 && nH <= max_size )
                if PI(inA,inC,type) ~= 0
                    PI_type_temp(type,nH) = PI_type_temp(type,nH) + PI(inA,inC,type);
                end
            end
        end
    end
end
% No need to normalise, as it is already a distribution
PI_type = PI_type_temp;
