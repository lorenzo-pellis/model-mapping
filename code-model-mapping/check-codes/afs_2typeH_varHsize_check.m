function afs_each_type = afs_2typeH_varHsize_check(H,PI,F,muH,NGM_G,Rh,shapeHvector,eta)

% This function is a variant of "afs_2typeH_varHsize" that uses older
% variants of the functions needed to compute the average final size of an
% epidemic with 2 types (adults and children) and a household composition
% distribution.
% 
% Difference from "afs_2typeH_varHsize":
%   - The self-consistent equations for the average final size try to mimic
%   what done in Ball, Britton and Sirl (2011), Equations (2.2) and (2.3),
%   rather than Equation (5.2) of Ball and Lyne (2001), but I don't really
%   understand the explanation of Ball, Britton and Sirl (they have a 
%   slightly different setup from mine), so I fiddled with the equations 
%   till I got the right results (it matches the output of the stochastic 
%   simulation)
%   - The equations of Ball, Britton and Sirl don't require the use of the
%   composition-biased household distribution, so I don't need to use PI as
%   an input here.
%   - The function I use here for computing the within-household epidemic
%   actually allows also initial infectives, which in my case are 0
% 
% Input:
%   - H: an (n+1)x(n+1) matrix with element (i,j) giving the probability 
%   that a randomly selected household has (i-1) adults and (j-1) children,
%   where n is the largest household size.
%   - F = [F_a,F_c], giving the fraction of adults and children in the
%   population (can be computed from H, but I've already done it before, so
%   I might as well pass it directly).
%   - muH: average size of a randomly selected household
%   - NGM_G: next-generation matrix of global contacts
%   - Rh: the mean number of infectious contact an infectives
%   makes with all other household members, throughout his/her entire
%   infectious period (the within-household R0)
%   - shapeHvector = [ h_ratio, thetaH, psiH, phiH ], containing relevant 
%   household parameters:
%       - h_ratio: c_a/c_c, where c_a and c_c are, respectively, the number
%       of contacts an adult and a child make per unit of time
%       - thetaH: household assortativity, i.e. fraction of contacts a
%       child made with other children. This is always set to NaN, as
%       random mixing is assumed in the household, so that the
%       assortativity is household-composition-dependent (it takes the
%       value of the fraction of (other) children in the household) and
%       therefore needs to be computed for each household composition
%       - psiH: relative susceptibility of children in the household (same
%       as global, by assumption that it is a property of the individual)
%       - phiH: relative infectivity of children in the household (same
%       as global, by assumption that it is a property of the individual)
%   - eta: exponent of frequency dependent transmission, set at 1 by
%   default and never changed
%
% Output: a 2x1 vector [z_a;z_c], where z_t is the fraction of individuals
% of type t (adult or children) ultimately infected. Note that they are not
% fractions of the full population, but of each group, i.e. 0 <= z_t <= 1
%
% Note 1: this function relies on function "afs_Addy_2type_const", which
% computes the average final size (afs) of either type (2 types only) of a 
% within-household epidemic, expressed as a fraction of the initial 
% susceptibles of each type, following the equations at page 314 of Ball, 
% Britton and Sirl.
%
% Note 2: Compared to the Ball, Britton and Sirl paper, I find it more 
% natural to work with the NGM of global contacts, and with the second 
% index referring to the infector and the first to the infectee. In the 
% paper they use the opposite convention (swapping rows and columns) and 
% use as input a matrix of 1-to-1 transmission rates. 
% 
% Note 3: The self-consistent relationships for the final size in Equations
% (2.2) and (2.3) come from a slightly different set-up compared to mine, 
% so I instead find it more intuitive to follow Equation (5.2) of Ball and
% Lyne (2001), i.e. using the composition-biased household composition
% distribution.
%
% Reference: Supplementary Methods, Section 1.2.7 of
% Pellis, L et al (2020), Nature Communications
% 
% Methodological references: 
% Ball, Britton and Sirl (2011), Journal of Mathematical Biology
% Addy, Longini and Haber (1991), Biometrics
% Ball and Lyne (2001), Advances in Applied Probability
%
% Author: Lorenzo Pellis
% Last update: 31-10-2019 

[ ntypes, ~ ] = size(NGM_G); % It's equal to 2 now

LambdaG = diag( F.^(-1) ) * NGM_G; 
% This is what they use is the paper for global contacts, i.e. it is still 
% a matrix of total infection rates (i.e. from 1 person towards everyone 
% else, or 1-to-all), but the NGM also contained information about the 
% fraction of individuals of each type on which infectious contacts are 
% distributed. To invert this, I divide each row by such fractions (hence 
% the diag(F.^(-1))). 
% 
% Note that LambdaG is transposed with respect to the paper. 

z0 = ones(ntypes,1); % strating point for the iterative method
displacement = @(s) ( s - afs_Addy_2type_compbiased(H,PI,Rh,shapeHvector,esc_prob_multitype(F,LambdaG,s),eta) );
displacement_check = @(s) ( muH*(s.*F') - afsn_Addy_2type_distr_initinf_check(H,Rh,shapeHvector,[0,0],esc_prob_multitype(F,LambdaG,s),eta) );
% Now I need to iteratively find the solution: in 1 dim I would use fzero,
% but in multiple dimensions I need to use fsolve
afs_each_type = fsolve( displacement, z0 );
afs_each_type_check = fsolve( displacement_check, z0 );

if sum(afs_each_type-afs_each_type_check)>1e-5
    disp(afs_each_type-afs_each_type_check);
    error('The two methods are giving different answers)');
end
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
%         afs_each_type = fsolve( displacement, z0, optimoptions(@fsolve,'Diagnostics','off','Display','off') );
%         afs_each_type = fsolve( displacement, z0, optimset('Diagnostics','off','Display','none') );
        afs_each_type = fsolve( displacement, z0 );
    end
end
