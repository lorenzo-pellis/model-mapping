function afsn = afsn_Addy_2type_const_initinf_check(h,i,Lh1to1,pi)

% This function is an older version of "afsn_Addy_2typeH_const". An even
% older version was computing the binomial coefficients when needed, rather
% than precomputing them, but that's way too inefficient to keep using that
% 
% Difference from "afsn_Addy_2typeH_const":
%   - the output is in numbers of individuals of each type, rather
%   than fractions
%   - The function I use here for computing the within-household epidemic
%   actually allows also initial infectives
% 
% This function computes the final size distribution of a within-household
% epidemic where individuals also escape global infection from outside.
% Here, it is restricted to 2 types (adults and children) and assumes a
% constant total individual infectivity (constant duration of infectious
% period, or a non-random time-since-infection infectivity profile).
% From the final outcome distribution, the average epidemic size for either
% type is calculated, as numbers of individuals of each type.
%
% Input:
%   - h = [h_a,h_c]: a vector with the household composition (h_a adults
%   and h_c children)
%   - i = [i_a,i_c]: a vector with the number of initial infectives of
%   either type
%   - Lh1to1: a matrix containing the 1-to-1 total infectivities, i.e.
%   element (i,j) gives the total infectivity a specified infective of type
%   j exerts on a specified individual of type i
%   - pesc = [q_a,q_c]: a vector with the probabilities of susceptible
%   adults and children, respectively, escaping global infection from 
%   outside the household.
% 
% Output: a 2x1 vector giving the average number of individuals of either  
% type t (adult or children) ultimately infected in the within-household  
% epidemic, i.e. they are values between 0 and h_t.
% 
% Note 1: The matrix Lh1to1 is transposed compared to the paper below,
% and its elements have been calculated outside the function, after 
% choosing the denominator n^eta or (n-1)^eta for the frequency dependent 
% transmission.
% Note 2: For efficiency, the binomial coefficients are pre-computed 
% outside this function, with nchoosek(i,j) stored in bincoeffmat(i+1,j+1)
% (to allow for 0 values for i and j).
% 
% Reference: Supplementary Methods, Section 1.2.7 of
% Pellis, L. et al (2020), Nature Communications
%
% Methodological references: 
% Ball, Britton and Sirl (2011), Journal of Mathematical Biology
% Addy, Longini and Haber (1991), Biometrics
%
% Author: Lorenzo Pellis
% Last update: 31-10-2019 

global bincoeffmat;

n = sum(h); % household size
s = h - i; % number of initial susceptibles
max = (s(1)+1)*(s(2)+1);
lm = zeros(1,max); % these are called hM in the paper
ls = zeros(1,max); % these are called hS in the paper
A = zeros(max);
b = zeros(max,1);

for v = 1:max
    v_pair = map_2typeH_1to2(s,v);
    av = v_pair(1);
    cv = v_pair(2);
    lm(v) = mgf( ( s(1) - av ) * Lh1to1(1,1) + ( s(2) - cv ) * Lh1to1(2,1), 0, 1, 1 );  % Lh is transposed with respect to the paper
    ls(v) = mgf( ( s(1) - av ) * Lh1to1(1,2) + ( s(2) - cv ) * Lh1to1(2,2), 0, 1, 1 );  % Lh is transposed with respect to the paper
    for w = 1:v
        w_pair = map_2typeH_1to2(s,w);
        aw = w_pair(1);
        cw = w_pair(2);
        if (aw<=av) && (cw<=cv)
            A(v,w) = bincoeffmat(s(1)-aw+1,av-aw+1) * bincoeffmat(s(2)-cw+1,cv-cw+1) /  ( lm(v)^(aw+i(1)) * ls(v)^(cw+i(2)) * pi(1)^(s(1)-av) * pi(2)^(s(2)-cv) );
        end
    end
    b(v) = bincoeffmat(s(1)+1,av+1) * bincoeffmat(s(2)+1,cv+1);
end
p = A \ b;

afs_temp = zeros(2,1);
for t = 1:max
    pair = map_2typeH_1to2(s,t);
    afs_temp = afs_temp + pair' * p(t); % It should be the same as above
end
afsn = afs_temp; % Average final size (afs) as numbers (so afsn), without counting the initial infectives
    
    
    
    
