function PI = create_2type_size_biased_distr(H)

% This function generates the composition-biased household composition 
% distribution, i.e. the composition distribution of the household a 
% randomly selected adult/child.
% 
% Input: an (n+1)x(n+1) matrix with element (i,j)
% giving the probability that a randomly selected household has (i-1)
% adults and (j-1) children, where n is the largest household size.
% 
% Output: a (n+1)x(n+1)x2 double matrix. The first layer has elements (i,j)
% giving the probability that the household of a randomly selected adult
% has (i-1) adults and (j-1) children. The second layer for the household
% of a randomly selected child. Note: row 1 in layer 1 is 0 and column 1 in
% layer 2 is 0.
% 
% Reference: Supplementary Tables 3 (3A for (:,:,1) and 3B for (:,:,2)) of
% Pellis, L. et al (2020), Nature Communications
% 
% Author: Lorenzo Pellis
% Last update: 29-12-2019

[ max_inA,max_inC ] = size(H); % inA = index of number of adults (starting from 1 for nA = 0), so inA = nA-1
PItemp = zeros(max_inA,max_inC,2); % 2 layers: first for the house of a randomly selected adults, second for a randomly selected child
PI = zeros(max_inA,max_inC,2);
part_sum = zeros(1,1,2);
for type = 1:2
    for inA = 1:max_inA
        nA = inA - 1;
        for inC = 1:max_inC
            nC = inC - 1;
            nH = [nA,nC];
            PItemp(inA,inC,type) = nH(type) * H(inA,inC);
            part_sum(type) = part_sum(type) + PItemp(inA,inC,type);
        end
    end
    PI(:,:,type) = PItemp(:,:,type)/part_sum(type);
end
