function  R0 = R0AH(chains,Rg)

% This function computes R0 for a model with households and 2 types (in my
% case adults and children). There is a variant of this function that also
% computes the fractions of adults and children from the dominant
% eigenvector and the NGM, but this is faster if you only need R0.
%
% Input:
%   - chains: this is a 2xN matrix of the average number of cases of either
%   type (one per row) in each generation. The average is over the
%   stochastic dynamics and the household size distribution, so N is the
%   largest possible number of generations of all household compositions
%   allowed by the population structure.
%   - Rg is the next generation matrix for global infections, i.e. Rg(1,2)
%   is the average number of type 1 (adults) infected by a single type 2
%   (child) through global contacts only, during the child's entire 
%   infectious period, in a fully susceptible population 
%
% Note 1: the Mmatrix is 2N x 2N and its particular structure is described
% in the main paper.
% 
% Note 2: this function relies of "get_dominant_eigenpair" to obtain the
% dominant eigenvalue and relative eigenvector if the big matrix
%
% Reference: Supplementary Methods, Section 1.2.4 of
% Pellis, L et al (2019), Nature Communications
%
% Author: Lorenzo Pellis
% Last update: 14/05/2019 

[ ntypes, N, ntypes_check ] = size(chains);
if ( ntypes ~= ntypes_check )
    error('Careful! The chain array does not have the right dimensions!')
end

alpha = zeros(N,2,2); % These all the non-zero coefficients appearing in Mmatrix
for tee = 1:2 % tee is short for infectee
    for tor = 1:2 % tor is short for infector
        for iN = 1:N
            alpha(iN,tee,tor) = Rg(tee,1)*chains(1,iN,tor) + Rg(tee,2)*chains(2,iN,tor);
        end
    end
end

Mmatrix = [      alpha(:,1,1)',            alpha(:,1,2)'          ;
            eye(N-1), zeros(N-1,1),       zeros(N-1,N)          ;
                 alpha(:,2,1)',            alpha(:,2,2)'          ;
                 zeros(N-1,N),         eye(N-1), zeros(N-1,1)   ];

R0 = max(eig(Mmatrix));
