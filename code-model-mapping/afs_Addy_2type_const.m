function afs = afs_Addy_2type_const(h,Lh1to1,pesc)

% This function computes the final size distribution of a within-household
% epidemic where individuals also escape global infection from outside.
% Here, it is restricted to 2 types (adults and children) and assumes a
% constant total individual infectivity (constant duration of infectious
% period, or a non-random time-since-infection infectivity profile).
% From the final outcome distribution, the average epidemic size for either
% type is calculated, as a fraction of the total individuals of each type.
%
% Input:
%   - h = [h_a,h_c]: a vector with the household composition (h_a adults
%   and h_c children)
%   - Lh1to1: a matrix containing the 1-to-1 total infectivities, i.e.
%   element (i,j) gives the total infectivity a specified infective of type
%   j exerts on a specified individual of type i
%   - pesc = [q_a,q_c]: a vector with the probabilities of susceptible
%   adults and children, respectively, escaping global infection from 
%   outside the household.
% 
% Output: a 2x1 vector giving the fraction of individuals of either type t 
% (adult or children) ultimately infected in the within-household epidemic. 
% The fractions are out of the total number of either type, i.e. they are
% values between 0 and 1.
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
% Ball and Lyne (2001), Advances in Applied Probability
%
% Author: Lorenzo Pellis
% Last update: 31-10-2019 

global bincoeffmat;

n = sum(h); % total household size
k = h(1); % number of adults
max = (h(1)+1)*(h(2)+1); % maximum index size for the household composition
lm = zeros(1,max); % These are called hM in the paper
ls = zeros(1,max); % These are called hS in the paper
% The variables above are where the assumptions on the duration of the
% infectious period comes in, through the moment-generating function
A = zeros(max); % Matrix of coefficients: A*p = b
b = zeros(max,1); % Vector on the right-hand side of the system to solve 

for v = 1:max % First I loop through the equation in the solution system
    v_pair = map_2typeH_1to2(h,v); % Pair of subscripts associated with index v
    av = v_pair(1); % Number of adults in composition indexed by v
    cv = v_pair(2); % Number of children in composition indexed by v
    % I'm using the following notation (household composition is (k,n-k))
    % to ensure I'm doing the calculations correctly as in the paper. Also,
    % notice the matrix is transposed compared to the paper. The other
    % parameters are just to use a constant infectious period with mean 1
    lm(v) = mgf( ( k - av ) * Lh1to1(1,1) + ( n - k - cv ) * Lh1to1(2,1), 0, 1, 1 );
    ls(v) = mgf( ( k - av ) * Lh1to1(1,2) + ( n - k - cv ) * Lh1to1(2,2), 0, 1, 1 );
    for w = 1:v % For each equation in the solution system I need to sum over all the possible final states
        w_pair = map_2typeH_1to2(h,w);
        aw = w_pair(1);
        cw = w_pair(2);
        if (aw<=av) && (cw<=cv)
            A(v,w) = bincoeffmat(k-aw+1,av-aw+1) * bincoeffmat(n-k-cw+1,cv-cw+1) /  ( lm(v)^aw * ls(v)^cw * pesc(1)^(k-av) * pesc(2)^(n-k-cv) );
        end
    end
    b(v) = bincoeffmat(k+1,av+1) * bincoeffmat(n-k+1,cv+1);
end
p = A \ b; % Probability of each possible final outcome

afs_temp = zeros(2,1); % Average number of cases infected of each type
for t = 1:max
    pair = map_2typeH_1to2(h,t);
    afs_temp = afs_temp + pair' * p(t);
end
assert(sum(afs_temp(h==0))==0); % Check that if there are no adults or no children in the household, the average number of cases of that type must be 0
afs = afs_temp ./ h'; % Average final size, as proportion of initial susceptibles of each type
afs(h==0) = 0; % If there are no adults or no children in the household, ensure we set the final size to 0 (rather than NaN from previous operation)
    
    
