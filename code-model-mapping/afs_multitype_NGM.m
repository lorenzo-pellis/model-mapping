function afs_each_type = afs_multitype_NGM(h,NGM)

% This function calculates the asymptotic average final size (afs), 
% expressed as a fraction of the total population, for a multitype model
% with an arbitrary number of types.
% 
% Inputs:
%   - h: a vector of proportions of different types in the population
%   (expected to sum up to 1, but if not, they are renormalised) 
%   - NGM: the next-generation matrix, where element (i,j) gives the total
%   number of cases of type i infected by a single case of type j,
%   throughout the entire infectious period, in a fully susceptible
%   population.
% 
% Output: a vector z with the fraction of individuals of each type who are
% ultimately infected.  Note that they are not fractions of the full 
% population, but of each group, i.e. 0 <= z_i <= 1
% 
% Reference: Supplementary Methods, Sections 1.3 of
% Pellis, L. et al (2020), Nature Communications
% 
% Author: Lorenzo Pellis
% Last update: 29-12-2019 
 
ntypes = length(h); % Number of types
if ~(sum(h)==1)
    warning('Warning: population group fractions don''t sum to 1, so I''m renormalising...')
    h = h / sum(h); % Renormalise into fractions of total population
end
[check1,check2] = size(NGM);
if check1 ~= check2
    error('Error: next-generation matrix non squared!')
end
if check1 ~= ntypes
    error('Error: number of types does not agree with the size of the next generation matrix!')
end

LambdaG = diag( h.^(-1) ) * NGM; 
% To mirror what done in the reference above, I compute the escaping 
% probabilities, from the matrix with 1-to-all total infectivities (i.e. 
% from 1 person towards everyone else), rather than the NGM, which also 
% contains information about the  fractions of individuals of each type on 
% which infectious contacts are distributed, i.e. Lambda = (l_{ij}) but 
% NGM = ( h_i l_{ij} ). To invert this, I divide each row by such fractions
% (hence the diag(F.^(-1))). 
% 
% Note that LambdaG is transposed with respect to the paper. 

z0 = ones(ntypes,1); % starting point for the iterative method
displacement = @(s) ( 1 - s - esc_prob_multitype(h,LambdaG,s) );
% Now I need to iteratively find the solution: in 1 dim I would use fzero,
% but in multiple dimensions I need to use fsolve
afs_each_type = fsolve( displacement, z0 );
% If the method has not converged or gives me negative solutions for some 
% reasons, I attempt to re-apply it again multiple times from different 
% starting points (halving both every time).
if ( ( afs_each_type(1) <= 0 ) || ( afs_each_type(2) <= 0 ) )
    max_rep = 10;
    n_rep = 1;
    while ( ( afs_each_type(1) <= 0 ) || ( afs_each_type(2) <= 0 ) )
        n_rep = n_rep + 1;
        if  ( n_rep > max_rep ) 
            warning([ 'Careful! No positive final size in afs_multitype_MRC even after ',num2str(max_rep),' initial starting points. I give the 0 solution' ]);
            afs_each_type = zeros(ntypes,1);
            break;
        end
        z0 = z0 / 2;
        afs_each_type = fsolve( displacement, z0 );
        % For debugging purposes:
        % afs_each_type = fsolve( displacement, z0, optimoptions(@fsolve,'Diagnostics','off','Display','off') );
        % afs_each_type = fsolve( displacement, z0, optimset('Diagnostics','off','Display','none') );
    end
end