function  [ R0, v_AH ] = R0AH_and_vAH(chains,Rg)

% This function computes R0 and the fractions of adults and children in 
% each generation for a model with households and 2 types (in my case 
% adults and children). This is the main function I use in my code, though 
% there are alterantive versions for cross-checking ("R0AH_vAH_and_NGM",
% which additionally attempts to compute an NGM - not sure it's right - and
% "R0AH_vAH_and_NGM_check", which computes all three outputs from a 
% different big matrix, constructed in a slower way; also "R0AH" only
% computes R0)
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
% Reference: Supplementary Methods, Sections 1.2.4 and 1.2.5 of
% Pellis, L. et al (2020), Nature Communications
%
% Author: Lorenzo Pellis
% Last update: 23-05-2019 

[ ntypes, N, ntypes_check ] = size(chains);
if ( ntypes ~= ntypes_check )
    disp('Careful! The chain tensor seem not have the right dimensions!');
end

alpha = zeros(N,2,2); % This gives all the coefficients appearing in the big matrix
for tee = 1:2
    for tor = 1:2
        for NN = 1:N
            alpha(NN,tee,tor) = Rg(tee,1)*chains(1,NN,tor) + Rg(tee,2)*chains(2,NN,tor);
        end
    end
end

Mmatrix = [      alpha(:,1,1)',            alpha(:,1,2)'          ;
            eye(N-1), zeros(N-1,1),       zeros(N-1,N)          ;
                 alpha(:,2,1)',            alpha(:,2,2)'          ;
                 zeros(N-1,N),         eye(N-1), zeros(N-1,1)   ];

% To reconstruct the proportions of adults and children in each generation,
% I need to reconstruct the number of adults and children in generation n
% from the eigenvector of Mmatrix, which contains information on primary
% cases of either type along generations n, n-1, n-2, etc. (rather than
% being all synchronous at generation n). To do so, I need to multiply them
% by the elements of the chains of within-household epidemic, and for each
% type of primary case, I can look at both adults and children in various
% generations, so v_short is 2*N, while v_long is 4*N
[ R0, v_short ] = get_dominant_eigenpair(Mmatrix);
T = [ diag([chains(1,:,1),chains(1,:,2)]); diag([chains(2,:,1),chains(2,:,2)]) ];
v_temp = T * v_short;
v_long = v_temp / sum(v_temp);
v_a = sum( v_long(1:2*N) );
v_c = sum( v_long((2*N+1):(4*N)) );
v_AH = [ v_a; v_c ];
