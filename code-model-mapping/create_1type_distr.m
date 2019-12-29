function H_single = create_1type_distr(H)

% This function generates the size distribution of a randomly selected 
% household.
% 
% Input: an (n+1)x(n+1) matrix with element (i,j)
% giving the probability that a randomly selected household has (i-1)
% adults and (j-1) children, where n is the largest household size.
% 
% Output: a 1xn vector, with i-th element giving the probability that 
% a randomly selected household has size i.
% 
% Reference: h_n in Supplementary Tables 2, 5 and 7 of 
% Pellis, L. et al (2019), Nature Communications
% 
% Author: Lorenzo Pellis
% Last update: 29-12-2019 

[ max_inA,max_inC ] = size(H);
if max_inA == max_inC
    H_single = zeros(1,max_inA-1);
else
    H_single = zeros(1,(max_inA-1)+(max_inC-1));
end
H_single_temp = H_single;
for inA = 1:max_inA
    for inC = 1:max_inC
        n = inA + inC - 2;
        if ( n~=0 && n<max_inA ) 
            H_single_temp(n) = H_single_temp(n) + H(inA,inC);
        end
    end
end
H_single = H_single_temp;
