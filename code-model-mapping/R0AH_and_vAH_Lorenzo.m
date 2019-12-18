function [ R0, v_AH, NGM ] = R0AH_vAH_and_NGM_check(chains,Rg)

% 31/01/2012

% Rg is the next generation matrix for global infections

% The Mmatrix is 4N x 4N and has been derived from the system of difference
% equations, but we need to be careful about dividing by 0. There are
% multiple possible choices. I choose one, switch to the other if the first
% has a division by 0. If the second has it too, it means that the
% numerator is 0 and I simply leave 0 in the Mmatrix.
% The result is correct, as it coincides with "R0AH_and_vAH", but
% it's more complicated.

[ ntypes, N, ntypes_check ] = size(chains);
if ( ntypes ~= ntypes_check )
    'Careful! The chain tensor does not have the right dimensions!'
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

