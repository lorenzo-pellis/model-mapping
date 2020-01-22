function z = afsH_varHsize_Britton_den(PI,Rg,Rh,distribution,iota,alpha,eta,den_str)

% This function calculates the average asymptotic final size (afs), 
% expressed as a fraction of the total population, for a single-type
% households model with a distribution of household sizes. Any distribution
% for the duration of the infectious period is allowed, as well as
% specification of how the 1-to-1 infectivity depends on the household
% size. This function relies on the 
% method described in the book by Andersson and Britton (2000), which is
% based on computing the within-household final size distribution with "a"
% initial infectives, which are those that escape infection from outside.
% This method is less efficient than finding directly the final size
% distribution with individuals escaping infection from outside and no
% initial infectives (using Addy, Longini and Haber (1991) or Ball and Lyne
% (2001)). Therefore I only use this function to cross-check the results
% more efficiently computed with the function "afsH_varHsize_Addy_den".
% 
% Input:
%   - PI: the size-biased household size distribution (single-type)
%   - Rg is the average number cases a single case makes through global 
%   contacts only, during the entire infectious period, in a fully 
%   susceptible population, i.e. the basic reproduction number for global 
%   contacts only
%   - Rh: the mean number of infectious contact an infectives
%   makes with all other household members, throughout his/her entire
%   infectious period (the within-household R0)
%   - distribution: this can be
%       - 0: for the constant iota (i.e. Dirac-delta-distributed r.v.)
%       - 1: for an exponential distribution, with mean iota
%       - 2: for a Gamma distribution, with mean iota and shape parameter alpha
%   - iota: mean of the r.v. (only for distributions 1 and 2)
%   - alpha: shape of the r.v. (only for distribution 2)
%   - eta: exponent of frequency dependent transmission, set at 1 by
%   default in the model mapping paper, but can be changed in general
%   - den_str: a string that specifies the form of frequency dependence:
%       - 'n': transmission has the form beta/n^eta
%       - 'n-1': transmission has the form beta/(n-1)^eta
% 
% Output: average final size (numbers)
% 
% Reference: Supplementary Methods, Section 1.4 of
% Pellis, L et al (2020), Nature Communications
%
% Methodological reference: book by Andersson & Britton (2000), Section 6.3
% 
% Author: Lorenzo Pellis
% Last update: 14-06-2019 

displacement = @(s) ( s - afsn_Britton_sizebiased_den(PI,exp(-Rg*s),Rh,distribution,iota,alpha,eta,den_str) );

% z = fzero( displacement, 1 ); % Start searching from 1, to get the largest solution (0 is always a solution)
% The line above is possibly inefficient. It's better to use a hand-made
% function "fzeroinfmax" that exploits the fact that 0 is an "inf" (it is
% always a solution but I don't want it) and 1 is a "max".
z = fzeroinfmax( displacement, 0, 1 ); 

