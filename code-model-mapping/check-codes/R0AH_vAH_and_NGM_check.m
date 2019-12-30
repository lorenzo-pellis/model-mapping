function [ R0, v_AH, NGM ] = R0AH_vAH_and_NGM_check(chains,Rg)

% This function computes R0, the fractions of adults and children in each 
% generation and the NGM corresponding to a pure 2-type model for a model 
% with households and 2 types (in my case adults and children). This is NOT
% the main function I use in my code, because I don't know whether the NGM 
% is correct: the main function is "R0AH_and_vAH". This is only used for
% cross-checking, as also is "R0AH_vAH_and_NGM", which uses a different
% big matrix and is faster. This function is slower because the matrix
% computed is 4Nx4N, rather than only 2Nx2N, and is constructed according
% to the intuitive approach in Section 3.2.2 of Pellis, Ball and Trapman
% (2012), extended to 2 types. I am not sure the extension is correct, but
% the result is the same
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
% Note 1: the Mmatrix is 4N x 4N and its particular structure is NOT 
% described in the main paper, but it's inspired by the reference below.
% 
% Note 2: this function relies of "get_dominant_eigenpair" to obtain the
% dominant eigenvalue and relative eigenvector if the big matrix
%
% Reference: Section 3.2.2 of Pellis, Ball and Trapman (2012), Mathematical
% Biosciences
% 
% Author: Lorenzo Pellis
% Last update: 23/05/2019 

[ ntypes, N, ntypes_check ] = size(chains);
if ( ntypes ~= ntypes_check )
    error('Careful! The chain tensor does not have the right dimensions!')
end

subMmatrix = zeros(N,N,2,2,2,2);
for tee = 1:2
    for tprimee = 1:2
        for tor = 1:2
            for tprimor = 1:2
                if tee == tprimee                    
                    for NN = 1:N
                        subMmatrix(1,NN,tee,tprimee,tor,tprimor) = Rg(tee,tor);
                    end
                end
            end
        end
        for NN = 2:N
            if chains(tprimee,NN-1,tprimee) == 0
                if chains(tee,NN,tprimee) ~= 0
                    subMmatrix(NN,NN-1,tee,tprimee,3-tprimee,tprimee) = chains(tee,NN,tprimee) / chains(3-tprimee,NN-1,tprimee);
                end
            else
                if chains(tee,NN,tprimee) ~= 0
                    subMmatrix(NN,NN-1,tee,tprimee,tprimee,tprimee) = chains(tee,NN,tprimee) / chains(tprimee,NN-1,tprimee);
                end
            end
        end        
    end
end

Mmatrix = [     subMmatrix(:,:,1,1,1,1),    subMmatrix(:,:,1,1,1,2),    subMmatrix(:,:,1,1,2,1),    subMmatrix(:,:,1,1,2,2)     ;
                subMmatrix(:,:,1,2,1,1),    subMmatrix(:,:,1,2,1,2),    subMmatrix(:,:,1,2,2,1),    subMmatrix(:,:,1,2,2,2)     ;
                subMmatrix(:,:,2,1,1,1),    subMmatrix(:,:,2,1,1,2),    subMmatrix(:,:,2,1,2,1),    subMmatrix(:,:,2,1,2,2)     ;
                subMmatrix(:,:,2,2,1,1),    subMmatrix(:,:,2,2,1,2),    subMmatrix(:,:,2,2,2,1),    subMmatrix(:,:,2,2,2,2)     ];
            
[ R0, v ] = get_dominant_eigenpair(Mmatrix);
v_a = sum( v(1:2*N) );
v_c = sum( v((2*N+1):(4*N)) );
v_AH = [ v_a; v_c ];

NGM = zeros(2);
for tor = 1:2
    v_in = zeros(size(v));
    v_temp = v_in;
    v_temp( ( (tor-1)*2*N + 1 ) : ( tor*2*N ) ) = v( ( (tor-1)*2*N + 1 ) : ( tor*2*N ) );
    if any(v_temp)
        v_in = v_temp / sum(v_temp);
    else
        v_in = v_temp;
    end
    v_temp = Mmatrix * v_in;
    for tee = 1:2
        v_out = zeros(size(v));
        for tprimor = 1:2
            v_out( ( (tee-1)*2*N + 1 ) : ( tee*2*N ) ) = v_temp( ( (tee-1)*2*N + 1 ) : ( tee*2*N ) );
        end
        NGM(tee,tor) = sum(v_out);
    end
end

