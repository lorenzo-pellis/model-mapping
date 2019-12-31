function [ R0, v ] = get_dominant_eigenpair(K)

% This function returns the dominant eigenvalue and the corresponding
% eigenvector of a primitive matrix. Perron-Frobenius theory guarantees
% that if a matrix with real entries is primitive, then it has only one
% eigenvalue with largest modulus, and such eigenvalue is real and
% positive. The corresponding eigenvector has all components of the same
% sign, so this code returns them all positive and normalised to sum to 1.
% A matrix with all positive entries is primitive. An irreducible matrix is
% not necessarily primitive, and can have 2 eigenvalues with the same
% modulus, so CHECK WHAT THE CODE DOES in this case (if of relevance), or
% in the case of reducible matrices.
% 
% Note: there is an alternative version of this function to extract only
% the dominant eigenvector ("get_dominant_eigenvector"). If interested only 
% in the dominant eigenvalue of A, just use max(eig(A)).
% 
% Author: Lorenzo Pellis
% Last update: 26/05/2019 

[ V, E ] = eig( K ); % First find all eigenvalues and eigenvectors
temp = max(E); % E is a matrix with eigenvalues on the diagonal, so max needs to be applied twice
[R0,index] = max(temp); % Finds the maximum and the index of it
v_temp = V(:,index); % Column vector
if ( all( v_temp>=-1e-10 ) || all( v_temp<=1e-10 ) ) % bypass numerical errors
    v = abs( v_temp / sum(v_temp) ); % the abs is not needed, except when there are numerical errors (leading to +0 and -0)
else
    warning('Problems: the eigenvector should have all positive elements! Return NaN');
    v = NaN(size(v_temp));
end
