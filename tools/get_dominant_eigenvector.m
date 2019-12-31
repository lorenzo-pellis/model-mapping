function v = get_dominant_eigenvector(K)

% This function returns the dominant eigenvector of a primitive matrix. 
% Perron-Frobenius theory guarantees that if a matrix with real entries is 
% primitive, then it has only one eigenvalue with largest modulus, and such
% eigenvalue is real and positive. The corresponding eigenvector has all 
% components of the same sign, so this code returns them all positive and 
% normalised to sum to 1. A matrix with all positive entries is primitive. 
% An irreducible matrix is not necessarily primitive, and can have 2 
% eigenvalues with the same modulus, so CHECK WHAT THE CODE DOES in this 
% case (if of relevance), or in the case of reducible matrices.
% 
% Note: there is an alternative version of this function to extract both
% the dominant eigenvalue and the corresponding eigenvector 
% ("get_dominant_eigenpair"). If interested only in the dominant eigenvalue
% of A, just use max(eig(A)).
% 
% Author: Lorenzo Pellis
% Last update: 26/05/2019 

[ V, E ] = eig( K ); % First find all eigenvalues and eigenvectors
temp = max(E); % E is a matrix with eigenvalues on the diagonal, so max needs to be applied twice
[~,index] = max(temp); % Finds the maximum and the index of it
v_temp = V(:,index); % Column vector
if ( all( v_temp>=-1e-10 ) || all( v_temp<=1e-10 ) ) % bypass numerical errors
    v = abs( v_temp / sum(v_temp) ); % the abs is not needed, except when there are numerical errors (leading to +0 and -0)
else
    warning('Problems: the eigenvector should have all positive elements!');
    v = v_temp; % Still return an output, to avoid getting stuck. This case will be rejected by other methods of control (I hope)
end
