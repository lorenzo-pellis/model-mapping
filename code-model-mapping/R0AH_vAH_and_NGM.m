function  [ R0, v_AH, NGM ] = R0AH_vAH_and_NGM(chains,Rg)

% This function computes R0, the fractions of adults and children in each 
% generation and the NGM corresponding to a pure 2-type model for a model 
% with households and 2 types (in my case adults and children). This is NOT
% the main function I use in my code, because I don't know whether the NGM 
% is correct: the main function is "R0AH_and_vAH". This is only used for
% cross-checking, as also is "R0AH_vAH_and_NGM_check", which uses a 
% different big matrix and is slower
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
% Pellis, L et al (2019), Nature Communications
%
% Author: Lorenzo Pellis
% Last update: 23/05/2019 

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

% CAREFUL! The rest of this function is an ATTEMPT to construct a 2x2 NGM 
% for adults and children, that somehow incorporates the complications of
% household structure. The method is untested and I'm not sure it does the 
% right thing, but the output is not used in the code.
% 
% The method spits out an NGM which has the same R0 and same v = [v_a;v_c]
% as it should, but I'm not sure I distribute the responsibility of
% infection correctly between the infectors.

resp_temp = zeros(N,2,2,2);
NGM = zeros(2);
for tprim = 1:2
    for tor = 1:2
        for tee = 1:2
            for NN = 1:(N-1)
                if chains(3-tprim,NN,tprim) ~= 0
                resp_temp(NN,tprim,tor,tee) = Rg(tee,tor) + chains(tee,NN+1,tprim) / sum(chains(:,NN,tprim));
                end
            end
            resp_temp(N,tprim,tor,tee) = Rg(tee,tor);
        end
    end
end
resp = reshape(resp_temp,[2*N,2,2]);
for tee = 1:2
    for tor = 1:2
        NGM(tee,tor) = resp(:,tor,tee)' * v_long(2*N*(tor-1) + (1:(2*N))) / sum(v_long(2*N*(tor-1) + (1:(2*N))));
    end
end
